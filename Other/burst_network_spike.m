

clear all; close all
% cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
% XLSdir contains excel files; needs to contain the list of recordings to
% include and this will be where the excel results files will be spat out
% .mat results files will be put in MATdir; MAT files containing spike
% matrices need to be in MATdir
sampling_fr = 25000;
code_dir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
% XLSdir = 'D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\Unsorted_byFileType\XLS_unsorted';
% MATdir = 'D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\Unsorted_byFileType\MAT_unsorted';
XLSdir = 'G:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\By_culture\2022_data';
MATdir = 'G:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\By_culture\2022_data';
% scriptsDir = 'D:\MECP2_2019_AD\Scripts_and_Output';
scriptsDir = 'C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
addpath(genpath(scriptsDir));
spike_suffix = '_cSpikes_L-0.0627_RP1.mat';
method = 'cwt';
threshold = 'L-0.0627';

% cd(XLSdir); inc_files = readtable('recs_included.xlsx');
% 
% % files = sortrows(files,'rec','ascend');
% 
% list = inc_files.file;
% for i = 1:length(list)
%     files(i,1).name = strcat(list{i} , spike_suffix);
% end
files = dir('M*Spikes*RP1.mat');

% find files with a certain phrase in the directory (redundant with above)
% files = dir('*MPT*cSpikes_L-0.0627_RF1.mat');

% Removes specific files
% files = files(~contains({files.name}, 'TTX','IgnoreCase',true));

samplingRate = sampling_fr; fs = sampling_fr;
% sync_win_s = 0.05;
%
% spikes_only = 0;
% burst_stats = 0;
% cor_ctrl = 0; % 1 means correlate randomly shuffled trains; be sure to change the sync window for the ctrl
% fc_figures = 0;%change to 1 to add plots
% g_figures = 0;%graph theory
% g_measures = 1;

%% within channel bursts
cd(MATdir)

warning('off','MATLAB:nearlySingularMatrix');
% bakkum method for individual channels
% N = 30; minChan = 3;
%                 [burstMatrix, burstTimes, burstChannels] = burstDetect(spikeMat, method, samplingRate,N,minChan);

fprintf(strcat('\n','\n','getting ISIn bursts within channels...','\n','\n'))
for i = 1:length(files)
    progressbar('elecs','files')
    filename = files(i).name
    load(filename);
    spikeMat = full(cSpikes);
    method = 'Bakkum'; % ask tim if bakkum method works with downsampling
    % in burstDetect there was a switch loop at the top that if using the
    % bakkum method sampling fr would be 25000
    minChan = 1;
    N = 10;
    active_chans = find(full(sum(cSpikes)) > N);
    for elec = 1 : size(spikeMat,2)
        %     spikeTrain = spikeMat(:,active_chans(elec));
        spikeTrain = spikeMat(:,elec);
        if sum(spikeTrain) >= N
            try
                [burstMatrix, burstTimes, burstChannels] = burstDetect(spikeTrain, method, samplingRate,N,minChan);
            catch
                burstMatrix     = 0;
                burstTimes      = 0;
                burstChannels   = 0;
            end
        else
            burstMatrix     = 0;
            burstTimes      = 0;
            burstChannels   = 0;
        end
        burstMatrices{elec} = burstMatrix;
        burstTimes_all{elec} = burstTimes;
        burstChannels_all{elec} = burstChannels;
        clear burstMatrix burstTimes burstChannels
        progressbar(elec/size(spikeMat,2),i/length(files))
    end
    progressbar(elec/size(spikeMat,2),i/length(files))
    
    
    % analyse bursts
    progressbar('elecs','files')
    for elec = 1:length(burstMatrices)
        %     numBursts(i,1) = length(burstMatrices{i});
        burstMatrix = burstMatrices{elec};
        burstTimes  = burstTimes_all{elec};
        if length(burstMatrix) > 0 & iscell(burstMatrix)
            for Bst=1:length(burstMatrix)
                sp_in_bst(Bst)=sum(sum(burstMatrix{Bst,1}));
                sp_times = find(burstMatrix{Bst,1}==1);
                sp_times2= sp_times(2:end);
                ISI_within = sp_times2 - sp_times(1:end-1);
                mean_ISI_w(Bst) = round(nanmean(ISI_within)/sampling_fr*1000,3); %in ms with 3 d.p.
                BLength(Bst) = size(burstMatrix{Bst,1},1)/sampling_fr*1000; %in ms
                within_bst_FR(Bst) = sp_in_bst(Bst) / BLength(Bst) *1000;
                clear ISI_within sp_times sp_times2 train
            end
            mean_num_sp_in_bst(elec)        = sum(sp_in_bst) / length(burstMatrix);
            total_num_sp_in_bst(elec)       = sum(sp_in_bst);
            burst_rate_elecs(elec)          = length(burstMatrix) / ((length(spikeMat)/fs)/60); %bursts/min
            mean_inBurst_FR(elec)           = nanmean(within_bst_FR);
            mean_ISI_w(elec)                = nanmean(mean_ISI_w);
            mean_BLength(elec)              = nanmean(BLength);
        else
            % disp('no bursts detected')
            mean_num_sp_in_bst(elec)    = NaN;
            total_num_sp_in_bst(elec)   = NaN;
            burst_rate_elecs(elec)      = NaN;
            mean_inBurst_FR(elec)       = NaN;
            mean_ISI_w(elec)            = NaN;
            mean_BLength(elec)          = NaN;
        end
        
        sp_times = find(spikeMat(:,elec)==1);
        sp_times2= sp_times(2:end);
        ISI_outside = sp_times2 - sp_times(1:end-1);
        
        mean_ISI_o(elec)        = round(nanmean(ISI_outside)/sampling_fr*1000,3); % in ms
        %         frac_spikes_inB(elec)   = total_num_sp_in_bst(elec) / sum(spikeMat(:,elec));
        progressbar(elec/length(burstMatrices),i/length(files))
        clear burstMatrix sp_in_bst
    end
    
    for elec = 1:length(burstMatrices)
        check_elec(1,elec) = ~isempty(burstMatrices{elec});
        check_elec(2,elec) = iscell(burstMatrices{elec});
    end
    bursting_electrodes = find(sum(check_elec) == 2);
    
    burstData(i).filename                = files(i).name;
    burstData(i).array_burstRate         = nanmedian(burst_rate_elecs(bursting_electrodes)); %bursts/min
    burstData(i).array_inBurstFR         = nanmedian(mean_inBurst_FR(bursting_electrodes)); %in Hz
    burstData(i).array_burstDur          = nanmedian(mean_BLength(bursting_electrodes)); %in ms
    %     burstData(i).array_fracInBursts      = nanmedian(frac_spikes_inB(bursting_electrodes));
    burstData(i).array_ISI_within        = nanmedian(mean_ISI_w(bursting_electrodes)); %in ms
    burstData(i).array_ISI_outside       = nanmedian(mean_ISI_o(bursting_electrodes)); %in ms
    burstData(i).array_fracInBursts      = (sum(total_num_sp_in_bst(~isnan(total_num_sp_in_bst)))  )  /  sum(sum(spikeMat));
    
    % get details of each burst for each rec. and elecrtode within that
    % rec.
    burstData(i).spike_matrices     = burstMatrices;
    burstData(i).burst_times        = burstTimes_all;
    %     below lines is not need because it is always one channel involved in
    %     the burst as it's within elec method!
    %     burstData(i,1).burst_channels     = burstChannels_all;
    
    %
    %                 %CV of IBI
    %                 %CV = st dev / mean
    %                 %get IBIs
    %                 end_times = burstTimes(1:end-1,2); %-1 because no IBI after end of last burst
    %                 sta_times = burstTimes(2:end,1); %start from burst start time 2
    %                 IBIs      = sta_times -end_times;
    %                 % calculate CV of IBI and non need to convert from samples to seconds
    %                 % (as relative measure it would be the same)
    %
    %                 % NOTE: these are based on the ISI across all channels!!!
    %                 output(file).mean_NBst_length_s = mean(NBLength);
    %                 output(file).num_Nbursts = length(burstTimes);
    %                 output(file).mean_chans_involved_in_Nbursts = mean(chans_involved);
    %                 output(file).mean_ISI_withinNbursts_ms  = mean(mean_ISI_w);
    %                 output(file).mean_ISI_outsideNbursts_ms = round(mean(ISI_outside)/sampling_fr*1000,3);
    %                 output(file).CVIofNBI = round((std(IBIs)/mean(IBIs)),3); %3 decimal places
    %                 % clear unneeded vars
    %                 clear end_times sta_times IBIs
    %                 output(file).NBurstRate=round(60*(nBursts/(length(spikeMat(:,1))/samplingRate)),3);
    %                 output(file).frac_in_Nburst = round(sp_in_bst/sum(sum(spikeMat)),3);
    %
    
    % burst rate /min
    
    % fraction spikes w/in bursts
    
    % average number of spikes within bursts
    
end


% log ISI
% spikeTimes = findSpikeTimes(downVec, 'seconds', samplingRate);

%% detect individual electrode bursts - Rank Surprise Method
% Gourévitch and Eggermont 2007
% not one of the four recommend by cotterill et al but performed not too
% badly in their results
progressbar('files')
clear burstMatrix electrodeBurst
fprintf(strcat('\n','\n','getting RS method bursts within channels...','\n','\n'))
for i = 1:length(files)
    filename = files(i).name;
    load(filename);
    spikeMat = full(cSpikes);
    downfactor = 25;
    downMatrix = downsample(spikeMat, downfactor);
    newfs = samplingRate/downfactor;
    spikeTimes = findSpikeTimes(spikeMat, 'seconds', newfs);
    
    numBurst = zeros(size(spikeMat, 2), 1);
    aveBurstSpikeNum = zeros(size(spikeMat, 2), 1);
    method = 'surprise';
    minChannel = 1; N = 3; %N here is min n spikes in a burst
    %     electrodeBurst = burstDetect(downMatrix, method, newfs,N,minChannel);
    %     %OUTPUT : archive_burst_RS - "rank surprise" values for each burst detected
    % %         archive_burst_length - burst length for each burst detected (in spikes)
    % %         archive_burst_start ? Spike number of burst start for each burst detected
    
    spikeTimes = findSpikeTimes(downMatrix, 'seconds', newfs);
    minSpike = N; % do not analyse trains with less than this many spikes
    for n = 1:length(spikeTimes)
        if length(spikeTimes{n}) < minSpike
            burstMatrix{n} = NaN;
        else
            [burstMatrix{n}.archive_burst_RS, ...
                burstMatrix{n}.archive_burst_length, ...
                burstMatrix{n}.archive_burst_start] = ...
                surpriseBurst(spikeTimes{n});
            if isempty(burstMatrix{n}.archive_burst_RS) % no bursts
                burstMatrix{n} = NaN;
            end
        end
    end
    electrodeBurst = burstMatrix;
    
    for n = 1:length(electrodeBurst)
        if ~ isfield(electrodeBurst{n}, 'archive_burst_RS')
            % check if there is a structure
            numBurst(n) = 0;
            aveBurstSpikeNum(n) = NaN;
        else
            numBurst(n) = length(electrodeBurst{n}.archive_burst_RS);
            % total number of burst from that electrode
            aveBurstSpikeNum(n) = nanmean(electrodeBurst{n}.archive_burst_length);
            % average number of spikes per burst from that electrode
        end
    end
    
    % burstTimes
    
    % total number of burst
    totalBurstCount(i) = sum(numBurst);
    arrayBurstCount(i) = nanmedian(numBurst);
    RS_burstData(i).totalBurstCount = totalBurstCount(i);
    RS_burstData(i).arrayBurstCount = arrayBurstCount(i);
    RS_burstData(i).totalBurstRate = totalBurstCount(i) / ((length(spikeMat)/samplingRate)/60) ;
    RS_burstData(i).arrayBurstRate = arrayBurstCount(i) / ((length(spikeMat)/samplingRate)/60) ;
    totalAveBurstSpikeNum(i) = nanmean(aveBurstSpikeNum);
    RS_burstData(i).totalAveBurstSpikeNum = totalAveBurstSpikeNum(i);
    RS_burstData(i).totalFracSpikesInBursts = (totalBurstCount(i) * totalAveBurstSpikeNum(i)) / sum(sum(spikeMat));
    % ^ estimate of total number of spikes in bursts across whole array
    RS_burstData(i).file = files(i).name;
    clear burstMatrix electrodeBurst numBurst aveBurstSpikeNum spikeTimes spikeMat downMatrix
    progressbar(i/length(files))
end
% Other statitics that may be useful:
% within burst firing frequency
% fraction of spikes w/in bursts
% CVofIBI
% number of electrodes activte during burst
% burst duration in ms, or frames

%% NNO burst detection method
% NNO is the name of the author: Nikolaas N. Oosterhof.
% The implementation is pulled from: https://github.com/nno/burstiDAtor.

%% firing regularity by
% Mochizuki et al.(2016)  source https://www.dropbox.com/s/7aiqgjgdymeow37/project_report_Tim_Sit.pdf?dl=0
progressbar('files')
fprintf(strcat('\n','\n','calculating firing regularity...','\n','\n'))
tic
for i = 1:length(files)
    filename = files(i).name;
    load(filename);
    spikeMat = full(cSpikes);
    spikeTimes = findSpikeTimes(spikeMat, 'seconds', samplingRate);
    spikeISI = findISI(spikeTimes);
    regularity = getReg(spikeISI, 'gamma', 10);
    totalReg(i) = nanmean(regularity);
    % freg = figure; histogram(log10(regularity));
    % freg2 = figure; scatter(numBursts, regularity);
    % freg3 = figure; r1 = regularity .* (1./numBursts);
    % r1(find(r1 == Inf)) = 0;
    % scatter(numBursts, r1);
    
    
    %{
shows very positively skewed
is the mean the best measure to take
there are no negative values but
the more positive the logshape,
the more regular the spiking activity,
= 0 reflects Poisson distributed spikes,
> 0 reflects regular firring, and
 < 0 reflects bursting activity
but bursting according to bakkum method correlates positively with reg.?
correlation is removed if you divide reg in each channel by numbursts
    %}
    FR_regularity(i).regularity = totalReg(i);
    FR_regularity(i).file = filename;
    progressbar(i/length(files))
    clear spikeMat regularity  spikeISI spikeTimes spikeMat filename
    toc
end


%% network spikes
folderPath ='D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output\';
networkSpikeDetectMethod = 'simple';
networkSpikeDetectParam.minChannel = 7;
networkSpikeDetectParam.binWidth = 0.1; % in seconds
fs = 25000;
jitterCriteria = true; % remove simultaneous spikes
clearvars -except files fs networkSpikeDetectMethod ...
    networkSpikeDetectParam jitterCriteria burstData ...
    RS_burstData FR_regularity method threshold
fprintf(strcat('\n','\n','identifying network spikes...','\n','\n'))
networkSpikeData = batchNetworkSpikeDetect(networkSpikeDetectMethod, networkSpikeDetectParam, fs, jitterCriteria,files);
%{
Here are what each value is:
        networkSpikeData{file, 1} = numNetworkSpike;
        networkSpikeData{file, 2} = networkSpikeTimes;
        networkSpikeData{file, 3} = participatingElectrodes;
        networkSpikeData{file, 4} = networkSpikeDetectParam.binWidth;
        networkSpikeData{file, 5} = networkSpikeDetectParam.minChannel;
        networkSpikeData{file, 6} = 'simple';
        networkSpikeData{file, 7} = fileNames(file);
%}
% saveFolder = folderPath;
% method = 'simple';
% batchVisualiseNetworkSpike(folderPath, networkSpikeData, method, saveFolder)
%current;y visualisation script doesnt work

%% save
disp('saving...')
fileName = strcat(method,'_' ,threshold, '_bursts_netwSpikes_regularity','','.mat');
% if loops below check if all four sections/measures have been calculated
% or only bakkum bursts
cd(MATdir)
if ~~exist('burstData') && ~~exist('RS_burstData') && ~~exist('FR_regularity') ...
        && ~~exist('networkSpikeData')
    save(fileName, 'burstData','RS_burstData','FR_regularity','networkSpikeData');
elseif ~~exist('burstData')
    fileName = strcat(method,'_' ,threshold, '_bursts','.mat');
    save(fileName, 'burstData');
elseif ~~exist('FR_regularity')
    fileName = strcat(method,'_' ,threshold, '_regularity','.mat');
    save(fileName, 'FR_regularity');
end

% NOTE: saving xls file in xls directory specified at the top
cd(XLSdir)
if ~~exist('burstData') && ~~exist('RS_burstData') && ~~exist('FR_regularity') ...
        && ~~exist('networkSpikeData')
    xldata1 = struct2table(burstData);
    xldata2 = struct2table(RS_burstData);
    xldata3 = struct2table(FR_regularity);
    xldata4 = cell2table(networkSpikeData);
    xldata4(:,'networkSpikeData2') = [];
    xldata4(:,'networkSpikeData3') = [];
    xldata4.Properties.VariableNames{2} = 'timeWindow';
    xldata4.Properties.VariableNames{1} = 'CountNetSpikes';
    xldata4.Properties.VariableNames{3} = 'MinNumChannels_in_a_netSpike';
    xldata4.Properties.VariableNames{4} = 'method';
    xldata4.Properties.VariableNames{5} = 'file';
    
    % xldata(:,'spikes') = [];
    % xldata(:,'thresholds') = [];
    % try
    %     xldata(:,'RC_node_IDs_allthresh') = [];
    %     xldata(:,'RC_node_IDs_maxthresh') = [];
    % catch
    % end
    writetable(xldata1(:,1:end-2), strcat(fileName(1:end-4),'BakkumBurstWithin'  ,'','.xlsx'))
    writetable(xldata2           , strcat(fileName(1:end-4),'RankSurpriseBursts' ,'','.xlsx'))
    writetable(xldata3           , strcat(fileName(1:end-4),'FiringRegularity'   ,'','.xlsx'))
    writetable(xldata4           , strcat(fileName(1:end-4),'NetworkSpikes'      ,'','.xlsx'))
elseif ~~exist('burstData')
    xldata1 = struct2table(burstData);
    writetable(xldata1(:,1:end-2), strcat(fileName(1:end-4),'BakkumBurstWithin','','.xlsx'))
elseif ~~exist('FR_regularity')
    xldata1 = struct2table(FR_regularity);
    writetable(xldata1(:,[2 1]), strcat(fileName(1:end-4),'FiringRegularity','','.xlsx'))
end

% save spike matrices for each burst for each electrode

