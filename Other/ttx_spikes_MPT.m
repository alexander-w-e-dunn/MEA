% get TTX spikes using absolute thresholds from corresponding electrodes

clear all; close all
rawDir = 'D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\By_culture\MPT\Recordings';
spikeDir = 'D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\By_culture\MPT\SpikesOnwards';
cd(rawDir)
% % organoid
% ttxfiles = dir('*FTD*TTX.mat');                   % where your .mat files are
% ttxfiles = ttxfiles(~contains({ttxfiles.name}, 'Spikes'));
% ttxfiles = ttxfiles(~contains({ttxfiles.name}, 'ttx4m'));
% ttxfiles = ttxfiles(~contains({ttxfiles.name}, '191210'));

% % MPT
ttxfiles = dir('*MPT*TTXAFTER*.mat');                   % where your .mat files are
ttxfiles = ttxfiles(~contains({ttxfiles.name}, 'Spikes'));

refPeriod_ms = 1;
method = 'abs';
L = 0;
original_ms = [3.5 4 4.5 5 5.5 6]; %multiplier used for mSpikes before TTX
% [log(0.001)/36.7368 log(0.01)/36.7368 log(0.1)/36.7368 log(10)/36.7368 log(100)/36.7368 log(1000)/36.7368]

progressbar('parameters','files','elecs')

for i = 1:length(original_ms)
    original_m = original_ms(i);
    for ttxfile = 1:length(ttxfiles)
        %strcat(ttxfiles(ttxfile).name(1:end-4), '_cSpikes_L',num2str(L), '.mat')
        ttxfileName = ttxfiles(ttxfile).name;
        ttxNameInd  = strfind(lower(ttxfileName),lower('TTX'));
        disp('getting thresholds from pre-ttx file')
        cd(spikeDir)
        try
            spikefile = strcat(ttxfileName(1:ttxNameInd-1),'TTXBEFORE','_mSpikes_',num2str(original_m),'_RF',num2str(refPeriod_ms),'.mat');
            load(spikefile,'thresholds');
        catch
            %manual correction for second ttx file
            if ~~contains(ttxfileName,'TTX2')
                spikefile = strcat(ttxfileName(1:end-10),'1_mSpikes_',num2str(original_m),'.mat');
                load(spikefile,'thresholds');
            else
                spikefile = strcat(ttxfileName(1:end-9),'1_mSpikes_',num2str(original_m),'.mat');
                load(spikefile,'thresholds');
            end
        end
        
        % loop to cancel process if file already done
        %     if ~~exist(strcat(ttxfiles(ttxfile).name(1:end-4), '_cSpikes_L',num2str(L), '.mat'))
        if ~isempty(dir(strcat(ttxfiles(ttxfile).name(1:end-4), '_aSpikes_based_on_mSpikes_',num2str(original_m), '*.mat')));
            disp('file done')
            
        else
            
            disp('loading ttx data')
            cd(rawDir)
            load(ttxfileName);
            
            pre_ttx_thresholds = thresholds;
            clear thresholds;
            
            disp('detecting spikes using absolute threshold')
            
            spikeMatrix = zeros(size(dat));
            finalData = zeros(size(dat)); %without initialising, you get memory error
            for elec = 1:length(pre_ttx_thresholds)
                
                multiplier = pre_ttx_thresholds(elec);
                
                % prevent spikes being detected in ref or grounded elecs
                max_spike_amplitude = -7;
                pre_ttx_thresholds(find(pre_ttx_thresholds > max_spike_amplitude)) = -100;
                
                [spikeMatrix(:, elec), finalData(:, elec), thresholds(elec)] = detectSpikes(dat(:, elec), method, multiplier,L,refPeriod_ms);
                
                %remove spikes from grounded electrodes
                spikeMatrix(:, find(pre_ttx_thresholds > max_spike_amplitude)) = 0;
                spikeMatrix(:, find(channels == 15)) = 0;
                progressbar(i/length(original_ms),ttxfile/length(ttxfiles),elec/length(pre_ttx_thresholds))
                
            end
            cd(spikeDir)
            % get spikes and save
            disp('saving ttx spike file')
            if strcmp(method,'cwt')
                cSpikes = sparse(spikeMatrix);
                %toc
                fileName = strcat(ttxfiles(ttxfile).name(1:end-4), '_cSpikes_L',num2str(L), '.mat');
                % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
                save(fileName, 'cSpikes','channels','thresholds','refPeriod_ms');
                
            elseif strcmp(method,'Manuel')
                mSpikes = sparse(spikeMatrix);
                %toc
                fileName = strcat(ttxfiles(ttxfile).name(1:end-4), '_mSpikes_',num2str(multiplier), '.mat');
                % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
                save(fileName, 'mSpikes','channels','thresholds','refPeriod_ms');
                
            elseif strcmp(method,'abs')
                aSpikes = sparse(spikeMatrix);
                %toc
                fileName = strcat(ttxfiles(ttxfile).name(1:end-4), '_aSpikes_based_on_mSpikes_',num2str(original_m),'_mthresh_',num2str(mean(pre_ttx_thresholds)), '.mat');
                save(fileName, 'aSpikes','channels','thresholds','refPeriod_ms');
            else
                disp('method inputted incorrectly')
            end
            
        end
        
        progressbar(i/length(original_ms),ttxfile/length(ttxfiles),0)
        
    end
    progressbar(i/length(original_ms),ttxfile/length(ttxfiles),0)
end