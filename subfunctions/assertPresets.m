function intPresetsCreated = assertPresets()
	% get target path
	intPresetsCreated = 0;
	cellExperiments = {'RunReceptiveFieldMapping','RunDriftingGratings','RunNaturalMovie','RunOptoStim'};
	
	%% generate data & pre-allocate
	for intExp=1:numel(cellExperiments)
		strExp = cellExperiments{intExp};
		try
			sStimPresets = loadStimPreset(1,strExp);
		catch
			if strcmp(strExp,'RunReceptiveFieldMapping')
				
				sStimPresets = struct;
				sStimPresets.strExpType = strExp;
				sStimPresets.dblSecsBlankAtStart = 3;
				sStimPresets.dblSecsBlankPre = 0.3;
				sStimPresets.dblSecsStimDur = 0.6;
				sStimPresets.dblSecsBlankPost = 0.1;
				sStimPresets.dblSecsBlankAtEnd = 3;
				
				intPresetIdx = saveStimPreset(sStimPresets,strExp);
				intPresetsCreated = intPresetsCreated + 1;
			elseif strcmp(strExp,'RunDriftingGratings')
				for intSet=1:3
					if intSet == 1
						sStimPresets = struct;
						sStimPresets.strExpType = strExp;
						sStimPresets.intNumRepeats = 100;
						sStimPresets.vecOrientations = 0:15:345;
						sStimPresets.vecOrientationNoise = 0;
						sStimPresets.dblSecsBlankAtStart = 3;
						sStimPresets.dblSecsBlankPre = 0.4000;
						sStimPresets.dblSecsStimDur = 1;
						sStimPresets.dblSecsBlankPost = 0.1000;
						sStimPresets.dblSecsBlankAtEnd = 3;
					elseif intSet == 2
						sStimPresets = struct;
						sStimPresets.strExpType = strExp;
						sStimPresets.intNumRepeats = 400;
						sStimPresets.vecOrientations = [0 5 90 95];
						sStimPresets.vecOrientationNoise = [0 2 0 2];
						sStimPresets.dblSecsBlankAtStart = 3;
						sStimPresets.dblSecsBlankPre = 0.4000;
						sStimPresets.dblSecsStimDur = 1;
						sStimPresets.dblSecsBlankPost = 0.1000;
						sStimPresets.dblSecsBlankAtEnd = 3;
						boolGenNoise = 1;
					elseif intSet == 3
						sStimPresets = struct;
						sStimPresets.strExpType = strExp;
						sStimPresets.intNumRepeats = 16;
						sStimPresets.vecOrientations = [0 45 90 135];
						sStimPresets.vecOrientationNoise = 0;
						sStimPresets.dblSecsBlankAtStart = 3;
						sStimPresets.dblSecsBlankPre = 0.4000;
						sStimPresets.dblSecsStimDur = 1;
						sStimPresets.dblSecsBlankPost = 0.1000;
						sStimPresets.dblSecsBlankAtEnd = 3;
						sStimPresets.vecStimulusSize_deg = [4 6 9 14 21 32 48 72];
					end
					intPresetIdx = saveStimPreset(sStimPresets,strExp);
					intPresetsCreated = intPresetsCreated + 1;
				end
			elseif strcmp(strExp,'RunNaturalMovie')
				sStimPresets = struct;
				sStimPresets.strExpType = strExp;
				sStimPresets.intNumRepeats = 200;
				sStimPresets.dblSecsBlankAtStart = 3;
				sStimPresets.dblSecsBlankPre = 0;
				sStimPresets.dblSecsBlankPost = 0;
				sStimPresets.dblSecsBlankAtEnd = 3;
				%sStimPresets.dblSecsStimDur = 500/60
				intPresetIdx = saveStimPreset(sStimPresets,strExp);
				intPresetsCreated = intPresetsCreated + 1;
			elseif strcmp(strExp,'RunOptoStim')
				%2) Pulse trains (as discussed)
				%--> 5 pulses, each 5ms long, varying intervals between the pulses according to the frequency
				%Frequencies = 1Hz, 2Hz, 5Hz, 10Hz, 20Hz, 50Hz, 100Hz
				%After each pulse train 10s pause
				%Can run up to 500 trials (as discussed it should repeat the frequencies in a row and then start at 1Hz again) or until stopped
				
				sStimPresets = struct;
				sStimPresets.strExpType = strExp;
				sStimPresets.dblPrePostWait = 2;%secs
				sStimPresets.intRepsPerPulse = 5;%count
				sStimPresets.intTrialNum = 40;%count
				sStimPresets.dblPulseWait = 2;%secs, at least ~0.2s
				sStimPresets.vecPulseITI = 1./[1 2 5 10 20 50];%secs
				dblPulseDur = 20/1000;%secs
				sStimPresets.vecPulseDur = dblPulseDur*ones(size(sStimPresets.vecPulseITI));%secs
				sStimPresets.dblPulseWaitSignal = sStimPresets.dblPulseWait/2;
				sStimPresets.dblPulseWaitPause = sStimPresets.dblPulseWait - sStimPresets.dblPulseWaitSignal;
				intPresetIdx = saveStimPreset(sStimPresets,strExp);
				intPresetsCreated = intPresetsCreated + 1;
			else
				error([mfilename ':InstallInconsistency'],sprintf('Preset installation error; "%s" is not recognized',strExp));
			end
		end
	end
end