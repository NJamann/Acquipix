function [vecSpikeCh,vecSpikeT,dblTotT] = DP_DetectSpikes(matData, sP, vecChanMap)
	%DP_DetectSpikes Performs fast GPU-accelerated thresholded-spike detection
	%	[vecSpikeCh,vecSpikeT,dblTotT] = DP_DetectSpikes(matData, sP, vecChanMap)
	%
	%input:
	% - matData [Ch x T]: gpuArray [Channel by Timepoint] data matrix
	% - sP [struct]: structure containing parameters, same as KiloSort2's ops
	%	.tstart:	time offset in bins to start spike detection (default: 0)
	%	.tend:		sample # to stop spike detection (default: inf)
	%	.fshigh:	high-pass cut-off (default: 150)
	%	.fs:		sampling frequency (default: 30000)
	%	.spkTh:		negative spike threshold (default: -6)
	%	.winEdge:	overlapping window edge (default: ceil((fs*2e-3)/2)*2+1)
	%	.NT:		batch size (default: 65600)
	%	.CAR:		switch for common average rereferencing; 0=off,1=mean,2=median (default: 1)
	%	.minType:	switch for running minimum; 1=fast,2=high-precision (default: 1)
	%	.fslow:		offset to start spike detection (optional, not recommended)
	% - vecChanMap [Ch x 1]: vector with which channels to use
	%
	%output:
	% - vecSpikeCh; uint16 gpuArray, channel origin of spike
	% - vecSpikeT; uint32 gpuArray, time of spike in ms after data start
	% - dblTotT; total processed time in seconds
	%
	%Based on Kilosort2's get_good_channels.m
	%
	%Version history:
	%1.0 - 5 Dec 2019
	%	Created by Jorrit Montijn
	
	%get parameters
	intStartT = getOr(sP,'tstart',0);
	intStopT = getOr(sP,'tend',inf);
	dblHighPassFreq = getOr(sP,'fshigh',150);
	dblSampFreq = getOr(sP,'fs',30000);
	dblSpkTh = getOr(sP,'spkTh',-6);
	intWinEdge = getOr(sP,'winEdge',ceil((dblSampFreq*2e-3)/2)*2+1);
	intBuffT = getOr(sP,'NT',65600);
	intCAR = getOr(sP,'CAR',1);
	intType = getOr(sP,'minType',1);
	if isfield(sP,'fslow') && ~isempty(sP.fslow) && sP.fslow<sP.fs/2
		dblLowPassFreq = sP.fslow;
	else
		dblLowPassFreq = [];
	end
	
	%build filter
	if ~isempty(dblLowPassFreq)
		[b, a] = butter(3, [dblHighPassFreq/dblSampFreq,dblLowPassFreq/dblSampFreq]*2, 'bandpass');
	else
		[b, a] = butter(3, dblHighPassFreq/dblSampFreq*2, 'high');
	end
	
	%get chan map
	if ~exist('vecChanMap','var') || isempty(vecChanMap)
		vecChanMap = 1:size(matData,1);
	end
	
	%define variables
	intNewFreq = 1000;
	dblFreqConverter = dblSampFreq/intNewFreq;
	vecSpikeCh = gpuArray.zeros(5e4,1, 'uint16');
	vecSpikeT = gpuArray.zeros(5e4,1, 'uint32'); %subsample to 1 kHz
	intSpikeCounter = 0;
	
	%check if gpuarray
	matData = gpuArray(matData);
	
	%get starting times
	vecStartBatches = (1 + intStartT):(intBuffT-2*intWinEdge):min((size(matData,2)-intBuffT),intStopT);
	
	% detect rough spike timings
	for intBatch=1:numel(vecStartBatches)
		intStart = vecStartBatches(intBatch);
		matBuffer = matData(vecChanMap,(intStart:(intStart+intBuffT-1)));
		
		% apply filters and median subtraction
		matFiltered = DP_gpufilter(matBuffer, b, a, intCAR);
	
		% very basic threshold crossings calculation
		matFiltered = matFiltered./std(matFiltered,1,1); % standardize each channel ( but don't whiten)
		matMins = DP_FindMins(matFiltered, 30, 1, intType); % get local minima as min value in +/- 30-sample range
		vecSpkInd = find(matFiltered<(matMins+1e-3) & matFiltered<dblSpkTh); % take local minima that cross the negative threshold
		[vecT, vecCh] = ind2sub(size(matFiltered), vecSpkInd); % back to two-dimensional indexing
		vecCh(vecT<intWinEdge | vecT>intBuffT-intWinEdge) = []; % filtering may create transients at beginning or end. Remove those.
		
		%save time of spike (xi) and channel (xj)
		if intSpikeCounter+numel(vecCh)>numel(vecSpikeCh)
			vecSpikeCh(2*numel(vecSpikeCh)) = 0; % if necessary, extend the variable which holds the spikes
			vecSpikeT(2*numel(vecSpikeCh)) = 0; % if necessary, extend the variable which holds the spikes
		end
		vecSpikeCh(intSpikeCounter + [1:numel(vecCh)]) = vecCh; % collect the channel identities for the detected spikes
		vecSpikeT(intSpikeCounter + [1:numel(vecT)]) = (vecT+intStart)/dblFreqConverter; % collect the channel identities for the detected spikes
		
		intSpikeCounter = intSpikeCounter + numel(vecCh);
	end
	
	%calculate outputs
	dblTotT = (intBuffT*intBatch)/dblSampFreq;
	vecSpikeCh = vecSpikeCh(1:intSpikeCounter);
	vecSpikeT = vecSpikeT(1:intSpikeCounter);
	%vecSpikeRatePerChannel = accumarray(vecSpikeCh,1) ./ dblTotT;
	