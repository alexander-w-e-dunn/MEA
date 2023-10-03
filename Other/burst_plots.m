% INPUTS
% spikes_fileName:  Single filename containing binary spikematrix
% bursts_fileName:  .mat file containing burst data for spikes_fileName
% spikesDir:        where the file, spikes_fileName.mat is saved
% burstsDir:             "          bursts_fileName      " 
% scriptsDir:       where this file is saved as well as other
%                   functions/scripts (.m files) used in this function
% sampling_rate:    samples/s (Hz) of recordings 
% 
% OUTPUTS
% figure showing burstMatrix (see function below for details)
% NOTE: currently the burstOnly_spikeMatrix does not work (it is just a
% matrix of zeros so can be deleted
% 
% example use:
% 
% clear all
% scriptsDir      = 'D:\MECP2_2019_AD\Scripts_and_Output';
% cd(scriptsDir)
% spikes_fileName = 'MPT190403_2B_DIV21_cSpikes_L-0.0627_RF1.mat';
% bursts_fileName = 'Bakkum_L-0.0627_bursts.mat';
% spikeDir        = 'D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\Unsorted_byFileType\MAT_unsorted';
% burstDir        = 'D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\Unsorted_byFileType\MAT_unsorted';
% sampling_rate   = 25000;
% [burstMatrix,downBstMat,burstOnly_spikeMatrix] = burst_plots(spikes_fileName, bursts_fileName,spikeDir,burstDir,scriptsDir,sampling_rate);

function [burstMatrix,downBstMat,burstOnly_spikeMatrix] = burst_plots(spikes_fileName, bursts_fileName,spikeDir,burstDir,scriptsDir,sampling_rate);

addpath(genpath(scriptsDir));
% get burst data and spike matrix
[all_burst_times,spikes,spikeTimes,sampling_rate]   = getBurstData(spikes_fileName, bursts_fileName,spikeDir,burstDir);
% get burst only spike matrix and burst matrix
[burstMatrix,burstOnly_spikeMatrix]                 = getBurstSpikeMatrices(spikes);
% burst-only spike matrix means all spikes outside of burst times are
% removed for each electrode
% burst matrix means for each electrode, there are 1s in every sample in
% which it is bursting and 0 when not bursting

% plot the first 1 s of the burstMatrix to show periods of bursting and not bursting
% f1 = figure; f1.Position = [700 413 520 372];
% c = gray;
% colormap(c([256 1], :));
% imagesc(burstMatrix(1:25000,:)');
% cb = colorbar; cb.Location = 'Southoutside';
% caxis( [0 1] ); % cb.Label.String = '';
% cb.Ticks = [0.25 0.75]; aesthetics; MEAgraphics(1);
% cb.TickLabels = {'Not-bursting','Bursting'};
% xlabel('Sample #'); ylabel('Electrode');

% plot the average number of channels bursting in consecutive 1 s time bins
% f2 = figure;
% plot the raster with histogram on top
f3 = figure; f3.Position = [600 42 712 954];
subplot(6,1,1)
% downsample to one ms
Hz                              = 1000; % downsample to 1 sample / ms
new_n_samples                   = (length(burstMatrix)/sampling_rate) * Hz;
% get mean of time bins
downBstMat                      = downSampleMean(burstMatrix, new_n_samples);
% set any bin containing a burst to 1
downBstMat(find(downBstMat>0))  = 1;
% get vector showing the number of channels in which bursting is occuring
% in each time bin
v                               = sum(downBstMat,2);
bin_duration_ms                 = 10000;
number_of_bins                  = length(downBstMat) / bin_duration_ms;
% sum bursting channels in each time bin
num_bars                        = length(v) / bin_duration_ms;
v_to_plot                       = downSampleMean(v,num_bars);
bar(1:length(v_to_plot) , v_to_plot ,'FaceColor','k','BarWidth',1) ; aesthetics;
ylabel('Num. electrodes'); xlabel('Time (s)');
% correct x axis
xt = xticks;
xticklabels(xt * (bin_duration_ms/1000) );
% turn off xlabels so they are only at bottom
xticklabels([]); xlabel([]);
% axis tight
hold on;
plot(1:length(v_to_plot),smooth(v_to_plot),'-r','LineWidth',1.5)

% plot burst raster underneath
subplot(6,1,2:6)
imagesc(downBstMat');
c = gray;
colormap(c([256 1], :));
cb = colorbar; cb.Location = 'Southoutside';
caxis( [0 1] ); % cb.Label.String = '';
cb.Ticks = [0.25 0.75]; aesthetics; MEAgraphics(1);
cb.TickLabels = {'Not-bursting','Bursting'};
xlabel('Time (s)'); ylabel('Electrode');
xt = xticks;
xticklabels(xt /1000);


% sparse burst matrix to save space
burstMatrix = sparse(burstMatrix);

%% function to get spike Matrix and burst matrix
    function [all_burst_times,spikes,spikeTimes,sampling_rate] = getBurstData(spikes_fileName, bursts_fileName,spikeDir,burstDir)
        
        % get spike matrix
        cd(spikeDir);       spikedata       = load(spikes_fileName);
        fields                              = fieldnames(spikedata);
        spikedata_cell                      = struct2cell(spikedata);
        % get number of samples/s and channel IDs
        try
            sampling_rate                       = spikedata.fs;
        catch
            sampling_rate                       = 25000;
            disp('WARNING: sampling rate assumed to be 25 kHz')
        end
        channels                            = spikedata.channels;
        % find spikes variable
        for f = 1:length(fields)
            if contains(fields{f},'spikes','IgnoreCase',true)
                spikedataindex = f;
            end
        end
        spikes = spikedata_cell{spikedataindex}; clear spikedata_cell
        
        % find burst data for this spike filename and load burst data
        cd(burstDir); load(bursts_fileName);
        for file = 1 : length(burstData)
            if burstData(file).filename == spikes_fileName
                burstdataindex = file;
            end
        end
        burstData   = burstData(burstdataindex);
        all_burst_times = burstData.burst_times;
        clear burstData % to save memory
        clear spikedata spikedata_cell
        
        % get spike times in samples for each channel
        for channel = 1 : size(spikes,2)
            spikeTimes{:,channel} = find(spikes(:,channel) == 1);
        end
    end

%% function to create a burst-only spike matrix and burst matrix
    function [burstMatrix,burstOnly_spikeMatrix] = getBurstSpikeMatrices(spikes)
        % burst-only spike matrix means all spikes outside of burst times are
        % removed for each electrode
        % burst matrix means for each electrode, there are 1s in every sample in
        % which it is bursting and 0 when not bursting
        
        % initialise BURST MATRIX
        burstMatrix = zeros(size(spikes));
        for electrode = 1 : size(spikes,2)
            bst_times = all_burst_times{1,electrode};
            if ~isempty(bst_times) && sum(bst_times(:)) ~= 0
                for burst_num = 1 : size(bst_times,1)
                    burstMatrix( bst_times(burst_num,1) : bst_times(burst_num,2) , electrode )   =   1;
                end
            end
        end
        
        % initialise BURST-ONLY SPIKE MATRIX
        burstOnly_spikeMatrix = spikes;
        for electrode = 1 : size(spikes,2)
            % get vector of spike times for electrode
            v_spikes = spikeTimes{electrode};
            % get vector of burst times for that electrode
            v_bursts = find(burstMatrix(:,electrode) == 1);
            % find and remove unique values (i.e. ones that are in
            % spike matrix but not in the burst matrix—i.e. spikes outside
            % of bursts
            spikeTimes_OutsideBursts = setdiff(v_spikes,v_bursts); % returns index of values in v_spikes but not in v_bursts
            % set times where there is a 1 in spikeMat and 0 in
            % burstMat—i.e spikes outside burst—and 0 in spikeMat and 1 in
            % burst—where there is no spike anyway—to 0.
            burstOnly_spikeMatrix(spikeTimes_OutsideBursts,electrode) = 0;
        end
    end
end