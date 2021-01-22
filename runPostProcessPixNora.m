% default options are in parenthesis after the comment
clear all;
%sites
cellRec{1}{1} = 'Z:\NIN202104_Long_range_connections\NIN202103_Neuropixels\Rawdata\SpikeGLX\20201208_NPX6_RunOptoNoraR01_g0';

matRunPrePro = [...
	1 1;...
    ];


%cellDepthCorrection{1}{1} = 0;%01

cellDepths{1}{1} = 0;

%cellDepths = cellfun(@plus,cellDepths,cellDepthCorrection);

cellMouseType{1}{1} = 'rbp4cre';

boolOnlyJson = false;

for intRunPrePro=[1]
	%% prepare
	% clear variables and select session to preprocess
	clearvars -except cellRec cellDepths cellMouseType matRunPrePro intRunPrePro boolOnlyJson
	vecRunPreProGLX = matRunPrePro(intRunPrePro,:);
	fprintf('\nStarting pre-processing of "%s" [%s]\n',cellRec{vecRunPreProGLX(1)}{vecRunPreProGLX(2)},getTime);
	
	%% set recording
	cellPath = strsplit(cellRec{vecRunPreProGLX(1)}{vecRunPreProGLX(2)},filesep);
	cellRecParts = strsplit(cellPath{end},'_');
	strMouse = cellRecParts{2};
	strExperiment = cellPath{end-1};
	strRecording = cellPath{end};
    cellRecSplit = strsplit(strRecording,'_');
	strExp = strcat('Exp',cellRecSplit{1}(1:4),'-',cellRecSplit{1}(5:6),'-',cellRecSplit{1}(7:8));
	strRecIdx = strcat('S',num2str(vecRunPreProGLX(1)),'L',num2str(vecRunPreProGLX(2))); %subject / location
	dblInvertLeads = true; %is ch1 deepest?
	dblCh1DepthFromPia = cellDepths{vecRunPreProGLX(1)}{vecRunPreProGLX(2)};
	strMouseType = cellMouseType{vecRunPreProGLX(1)}{vecRunPreProGLX(2)};
	
	% search for file
	strThisMouseIdx = getFlankedBy(strRecording,'_M','_','last');
	strThisRecIdx = getFlankedBy(strRecording,'R','_','last');
	if intRunPrePro == 4
		strThisRecIdx = getFlankedBy(strRecording,'_R','_','first');
		intLastPupilStop = 5;
	else
		intLastPupilStop = 1;
	end
	if vecRunPreProGLX(1) < 2
		strRec = ['*' strThisMouseIdx '*R' strThisRecIdx '*'];
	else
		strRec = ['*R' strThisRecIdx '*' strThisMouseIdx '*'];
	end

	%% set & generate paths
	strThisPath = mfilename('fullpath');
	strThisPath = strThisPath(1:(end-numel(mfilename)));
	strDataPath = strjoin(cellPath(1:(end-1)),filesep);
	strPathDataTarget = ['Z:\Jamann\PROJECTS - data\Neuropixels IVD202103\Preprocessed' filesep strExperiment filesep];
	if ~exist(strPathDataTarget,'dir'),mkdir(strPathDataTarget);end
	strChanMapFile = strcat(strThisPath,'subfunctionsPP\neuropixPhase3B2_kilosortChanMap.mat');
	strPathStimLogs = strcat(strjoin(cellPath(1:(end-2)),filesep),filesep,'Matlab');
	strPathEphys = fullfile(strDataPath,strRecording);
	fprintf('Processing recording at %s [%s]\n',strPathEphys,getTime);
	
	
	%% load NI sync stream times
	strFileNI = strcat(strRecording,'_t0.nidq.bin');
	fprintf('Loading syncing data %s [%s]\n',strFileNI,getTime);
	% Parse the corresponding metafile
	sMetaNI = DP_ReadMeta(strFileNI, strPathEphys);
	dblSampRateReportedNI = DP_SampRate(sMetaNI);
	intFirstSample = str2double(sMetaNI.firstSample);
	
	% Get NI data
	matDataNI = -DP_ReadBin(-inf, inf, sMetaNI, strFileNI, strPathEphys);
	[boolVecSyncPulses,dblCritValSP] = DP_GetUpDown(matDataNI(2,:));
	clear matDataNI;
	
	vecChangeSyncPulses = diff(boolVecSyncPulses);
	vecSyncPulseOn = (find(vecChangeSyncPulses == 1)+1);
	vecSyncPulseOff = (find(vecChangeSyncPulses == -1)+1);
	clear vecChangeSyncPulses boolVecSyncPulses;
	dblSampRateNI = mean(diff(vecSyncPulseOn));
	dblSampRateFault = (1-(dblSampRateReportedNI/dblSampRateNI));
	if dblSampRateFault < -1e-5 || dblSampRateFault > 1e-5
		error([mfilename 'E:SampRateFault'],sprintf('Sampling rate fault is high: %e. Please check!',dblSampRateFault));
	end
	
	%% load stimulus info
	%load logging file
	fprintf('Synchronizing multi-stream data...\n');
	dblLastStop = 0;
    strPathStimLogs = strcat(strPathStimLogs,filesep,strExp,filesep,strMouse);
	sFiles = dir(fullfile(strPathStimLogs,['*_' strMouse '*' strThisRecIdx '_*.mat']));
	intLogs = numel(sFiles);
	if intLogs == 0
		error([mfilename ':NoLogsFound'],'No log files found');
	else
		fprintf('\b   Found %d log files [%s]\n',intLogs,getTime);
	end
	
	%% determine temporal order
	cellFiles = {sFiles(:).name};
	vecTimes = nan(1,intLogs);
	for intLogFile = 1:intLogs
		cellSplit = strsplit(cellFiles{intLogFile}(1:(end-4)),'_');
		vecTimes(intLogFile) = str2double(cat(2,cellSplit{end-2:end}));
	end
	[dummy,vecReorderStimFiles] = sort(vecTimes);
	
	%% run
	cellStim = cell(1,intLogs);
	for intLogFile = 1:intLogs
		%% calculate stimulus times
		fprintf('>Log file "%s" [%s]\n',sFiles(vecReorderStimFiles(intLogFile)).name,getTime)
		cellStim{intLogFile} = load(fullfile(strPathStimLogs,sFiles(vecReorderStimFiles(intLogFile)).name));
		strStimType = cellStim{intLogFile}.structEP.strFile;
		
		intThisNumTrials = numel(~isnan(cellStim{intLogFile}.structEP.vecStimOffNI));
        cellStim{intLogFile}.structEP.ActOnNI = cellStim{intLogFile}.structEP.vecStimOnNI;
        cellStim{intLogFile}.structEP.ActOffNI = cellStim{intLogFile}.structEP.vecStimOffNI;
		if isfield(cellStim{intLogFile}.structEP,'ActOnNI') && ~all(isnan(cellStim{intLogFile}.structEP.ActOnNI))
			vecStimActOnNI = cellStim{intLogFile}.structEP.ActOnNI - intFirstSample/dblSampRateNI;
			vecStimActOffNI = cellStim{intLogFile}.structEP.ActOffNI - intFirstSample/dblSampRateNI;
		else
			%approximate timings
			vecStimOn = vecStimOnScreenPD/dblSampRateNI;
			vecStimOff = vecStimOffScreenPD/dblSampRateNI;
			%remove prior entries
			vecStimOn(vecStimOn < dblLastStop) = [];
			vecStimOff(vecStimOff < dblLastStop) = [];
			%ensure identical size
			if vecStimOff(1) < vecStimOn(1),vecStimOff(1) = [];end
			if numel(vecStimOn) > numel(vecStimOff),vecStimOn(end) = [];end
			%calc dur
			vecStimDur = vecStimOff - vecStimOn;
			
			%remove all durations shorter than single frame
			vecRemStims = vecStimDur <= cellStim{intLogFile}.structEP.dblStimFrameDur;
			vecStimDur(vecRemStims) = [];
			vecStimOn(vecRemStims) = [];
			vecStimOff(vecRemStims) = [];
			
			%get real but inaccurate timings
			vecStimActOnSecs = cellStim{intLogFile}.structEP.ActOnSecs - cellStim{intLogFile}.structEP.ActOnSecs(1);
			vecStimActOffSecs = cellStim{intLogFile}.structEP.ActOffSecs - cellStim{intLogFile}.structEP.ActOnSecs(1);
			
			%go through onsets to check which one aligns with timings
			vecSignalOnT = vecStimOnScreenPD/dblSampRateNI;
			intStims = numel(vecSignalOnT);
			vecError = nan(1,intStims);
			for intStartStim=1:intStims
				%select onsets
				vecUseSignalOnT = vecSignalOnT(intStartStim:end);
				
				%get ON times
				[vecStimOnTime,vecDiffOnT] = OT_refineT(vecStimActOnSecs,vecUseSignalOnT - vecUseSignalOnT(1),inf);
				vecError(intStartStim) = nansum(vecDiffOnT.^2);
			end
			[dblMin,intStartStim] = min(vecError);
			dblStartT = vecSignalOnT(intStartStim);
			
			%get probability
			vecSoftmin = softmax(-vecError);
			[vecP,vecI]=findmax(vecSoftmin,10);
			dblAlignmentCertainty = vecP(1)/sum(vecP);
			fprintf('Aligned onsets with %.3f%% certainty; start stim is at t=%.3fs\n',dblAlignmentCertainty*100,dblStartT);
			if (dblAlignmentCertainty < 0.9 || isnan(dblAlignmentCertainty)) && ~(intLogFile == 1 && intRunPrePro == 7) && ~(intLogFile == 2 && intRunPrePro == 10)
				error([mfilename 'E:CheckAlignment'],'Alignment certainty is under 90%%, please check manually');
			end
			%ensure same starting time
			vecStimActOnNI = vecStimActOnSecs + dblStartT;
			vecStimActOffNI = vecStimActOffSecs + dblStartT;
			
		end
		%remove missing stimuli
		vecRem = isnan(vecStimActOnNI) | isnan(vecStimActOffNI);
		dblLastStop = max(vecStimActOffNI);
		cellStim{intLogFile}.structEP = remStimAP(cellStim{intLogFile}.structEP,vecRem);
			
		
		% save to cell array
		cellStim{intLogFile}.structEP.ActOnNI = vecStimActOnNI;
		cellStim{intLogFile}.structEP.ActOffNI = vecStimActOffNI;
		cellStim{intLogFile}.structEP.SampRateNI = dblSampRateNI;
		
	end
	
	%% load clustered data into matlab using https://github.com/cortex-lab/spikes
	%load rez
	fprintf('Loading clustered spiking data at %s [%s]\n',strPathEphys,getTime);
	sLoad = load(fullfile(strPathEphys,'rez2.mat'));
	sRez = sLoad.rez;
	vecKilosortContamination = sRez.est_contam_rate;
	vecKilosortGood = sRez.good;
	
	% load some of the useful pieces of information from the kilosort and manual sorting results into a struct
	sSpikes = loadKSdir(strPathEphys);
	vecAllSpikeTimes = sSpikes.st;
	vecAllSpikeClust = sSpikes.clu;
	vecClusters = unique(vecAllSpikeClust);
	
	%get channel depth from pia
	sChanMap=load(strChanMapFile);
	vecChannelDepth = sChanMap.ycoords;
	vecChannelDepth = vecChannelDepth - max(vecChannelDepth);
	if dblInvertLeads,vecChannelDepth = vecChannelDepth(end:-1:1);end
	vecChannelDepth = vecChannelDepth + dblCh1DepthFromPia;
	
	%get cluster data
	fprintf('Assigning spikes to clusters... [%s]\n',getTime);
	[spikeAmps, vecAllSpikeDepth] = templatePositionsAmplitudes(sSpikes.temps, sSpikes.winv, sSpikes.ycoords, sSpikes.spikeTemplates, sSpikes.tempScalingAmps);
	vecAllSpikeDepth = dblCh1DepthFromPia - vecAllSpikeDepth;
	
	%remove nans
	for intStim=1:numel(cellStim)
		matStimOnOff = [cellStim{intStim}.structEP.vecStimOnNI;cellStim{intStim}.structEP.vecStimOffNI]';
		%remove nans
		vecRem = any(isnan(matStimOnOff),2);
		matStimOnOff(vecRem,:) = [];
		cellStim{intStim}.structEP = remStimAP(cellStim{intStim}.structEP,vecRem);
	end
	
	%% prepare spiking cell array
	intClustNum = numel(vecClusters);
	cellSpikes = cell(1,intClustNum);
	vecDepth = nan(1,intClustNum);
	for intCluster=1:intClustNum
		intClustIdx = vecClusters(intCluster);
		cellSpikes{intCluster} = vecAllSpikeTimes(vecAllSpikeClust==intClustIdx);
		vecDepth(intCluster) = mean(vecAllSpikeDepth(vecAllSpikeClust==intClustIdx));
	end
	if ~boolOnlyJson
	%% go through clusters
	sCluster = struct;
	for intCluster=1:intClustNum
		%get cluster idx
		intClustIdx = vecClusters(intCluster);
		vecSpikeTimes = cellSpikes{intCluster};
		sOut = getClusterQuality(vecSpikeTimes,0);
		
		%assign to object
		sCluster(intCluster).Exp = strExperiment;
		sCluster(intCluster).Rec = strRecording;
		sCluster(intCluster).Area = [];
		sCluster(intCluster).MouseType = strMouseType;
		sCluster(intCluster).Mouse = strMouse;
		sCluster(intCluster).Date = getDate;
		sCluster(intCluster).Depth = vecDepth(intCluster);
		sCluster(intCluster).Cluster = intCluster;
		sCluster(intCluster).IdxClust = intClustIdx;
		sCluster(intCluster).SpikeTimes = vecSpikeTimes;
		sCluster(intCluster).NonStationarity = sOut.dblNonstationarityIndex;
		sCluster(intCluster).Violations1ms = sOut.dblViolIdx1ms;
		sCluster(intCluster).Violations2ms = sOut.dblViolIdx2ms;
		sCluster(intCluster).Contamination = vecKilosortContamination(intCluster);
		sCluster(intCluster).KilosortGood = vecKilosortGood(intCluster);
		
		%msg
		fprintf('Cell %d/%d, Non-stat=%.3f, Viol=%.3f, Contam=%.3f [%s]\n',...
			intCluster,intClustNum,sOut.dblNonstationarityIndex,sOut.dblViolIdx2ms,vecKilosortContamination(intCluster),getTime);
	end
	
	%% load LFP data
	%{
	strFileLFP = strcat(strRecording,'_t0.imec0.lf.bin');
	fprintf('Filtering LFP data at %s [%s]\n',strFileLFP,getTime);
	sMetaLFP = DP_ReadMeta(strFileLFP, strPathEphys);
	matLFP = DP_ReadBin(0, inf, sMetaLFP, strFileLFP, strPathEphys, 'single');
	
	dblSampRateLFP = DP_SampRate(sMetaLFP);
	vecTimestampsLFP = (1:size(matLFP,2))/dblSampRateLFP;
	
	%filter each channel
	for intCh=1:size(matLFP,1)
		%get data
		vecFiltered = double(matLFP(intCh,:));
		
		%filter 50Hz
		vecWindow = [49.95 50.05]./(dblSampRateLFP./2);
		[fb,fa] = butter(2,vecWindow,'stop');
		vecFiltered = filtfilt(fb,fa,vecFiltered);
		
		%filter to 0.01-300Hz
		vecWindow2 = [0.01 300]./(dblSampRateLFP./2);
		[fb,fa] = butter(2,vecWindow2,'bandpass');
		vecFiltered = filtfilt(fb,fa,vecFiltered);
		matLFP(intCh,:) = cast(vecFiltered,'like',matLFP);
		
		%calc power
		%[vecFreq,vecPower] = getPowerSpectrum(vecFiltered,dblSampRateLFP,2);
		%loglog(vecFreq(5:end-4),conv(vecPower,normpdf(-4:4,0,2),'valid'));
	end
	%}
	end
	%% combine all data and save to post-processing data file
	%build Acquipix post-processing structure
	fprintf('Combining data and saving to disk... [%s]\n',getTime);
	sAP = struct;
	strFileOut = strcat(strExperiment,'_',strMouse,'_',strRecIdx,'_AP');
	strFileAP = strcat(strPathDataTarget,strFileOut,'.mat');
	strSecondPathAP = 'C:\DataAP\';
	strFileAP2 = strcat(strSecondPathAP,strFileOut,'.mat');
	%save LFP separately because of large size
	%sAP_LFP = struct;
	%strFileOutLFP = strcat(strFileOut,'_LFP');
	%strFileLFP = strcat(strPathDataTarget,strFileOutLFP,'.mat');
	
	%LFP
	%sAP_LFP.vecTimestampsLFP = vecTimestampsLFP;
	%sAP_LFP.matLFP = matLFP;
	%sAP_LFP.sMetaLFP = sMetaLFP;
	
	%stimulation & eye-tracking timings
	sAP.cellStim = cellStim;
	
	%probe data
	sAP.vecChannelDepth = vecChannelDepth;
	if ~boolOnlyJson
	%clusters & spikes
	sAP.sCluster = sCluster;
	
	%NI meta file
	sAP.sMetaNI = sMetaNI;
	%sAP.strFileLFP = strFileLFP;
	
	%save AP
	fprintf('Saving AP data to %s [%s]\n',strFileAP,getTime);
	save(strFileAP,'sAP');
	fprintf('Saving AP data to %s [%s]\n',strFileAP2,getTime);
	save(strFileAP2,'sAP');
	
	%save LFP
	%fprintf('Saving LFP data to %s [%s]\n',strFileLFP,getTime);
	%save(strFileLFP,'sAP_LFP','-v7.3');
	%fprintf('Done\n');!
	end
	
	%% generate json file for library
	%split recording name & define data
	cellData = strsplit(strRecording,'_');
	strRecDate = cellData{1};
	if ~exist('strFileLFP','var'),strFileLFP='';end
	
	%required fields
	sJson = struct;
	sJson.date = strRecDate;
	sJson.version = '1.0';
	sJson.project = 'MontijnNPX2020';
	sJson.dataset = 'Neuropixels data';
	sJson.subject = strMouse;
	sJson.investigator = 'Jorrit_Montijn';
	sJson.setup = 'Neuropixels';
	sJson.stimulus = 'VisStimAcquipix';
	sJson.condition = 'none';
	sJson.id = strjoin({strRecIdx,strMouse,strExperiment,strThisRecIdx},'_');
	
	%additional fields
	sJson.experiment = strExperiment;
	sJson.recording = strRecording;
	sJson.recidx = strRecIdx;
	sJson.mousetype = strMouseType;
	sJson.nstims = num2str(numel(cellStim));
	sJson.stims = strjoin(cellfun(@(x) x.structEP.strFile,cellStim,'uniformoutput',false),';');
	sJson.trials = strjoin(cellfun(@(x) num2str(numel(x.structEP.vecStimOnTime)),cellStim,'uniformoutput',false),';');
	sJson.nclust = numel(vecKilosortGood);
	sJson.ngood = sum(vecKilosortGood);
	
	%check meta data
	cellFields = fieldnames(sMetaNI);
	intMetaField = find(contains(cellFields,'recording'));
	if numel(intMetaField) == 1
		sJson.recording = strjoin({sMetaNI.(cellFields{intMetaField}),cellFields{intMetaField}},'_');
	else
		sJson.recording = '';
		warning([mfilename 'W:NoMetaField'],'Meta field not found in NI header file');
	end
	
	%file locations
	sJson.file_ap = strFileAP;
	sJson.file_ap2 = strFileAP2;
	sJson.file_lfp = strFileLFP;
	sJson.file_ni = sMetaNI.fileName;
	
	%save json file
	strJsonData = jsonencode(sJson);
	strJsonFileOut = strcat(strExperiment,'_',strMouse,'_',strRecIdx,'.json');
	strJsonTarget = fullfile(strPathDataTarget,strJsonFileOut);
	fprintf('Saving json metadata to %s [%s]\n',strJsonTarget,getTime);
	save(strJsonTarget,'strJsonData');
end