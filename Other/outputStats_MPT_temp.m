% check spikes for given set of recs ;
% input = .mat data files

%cd to data files folder
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.1.MPT190515_mat_files'

%get spike matrices (if already done comment out)
%batchGetSpike
clear all; close all
% gtdata = readtable(['C:\Users\alexd\OneDrive - University of Cambridge',...
%     '\Cam\PhD\Manuscripts\Mecp2\2022_data\cwt_L-0.0627_Mecp2_spikes_Nbursts_STTC_data.csv']);
scriptsDir      = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
addpath(genpath(scriptsDir));
%if got spikes already go to spikes folder:
% cd 'C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Manuscripts\Mecp2\2022_data'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\AnalysisOutput\By_culture\MPT\Spikes\Nano
% cd('body\mSpikes_5_RP1')
cd('C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\AnalysisOutput\By_culture\Collaboration_with_Paul')
% files = dir('*MPT*cSpikes_L-0.0627_RF1.mat');
files = dir('WT*Spikes*.mat');

files = files(~contains({files.name}, 'TTX','IgnoreCase',true));
files = files(~contains({files.name}, 'mspikes_6','IgnoreCase',true));
% exclude files that don't go up to div35
include = zeros(length(files),1);
clear related_files
for file = 1 : length(files)
    related_files = files(~~contains({files.name}, files(file).name(1 : strfind(files(file).name,'_DIV'))  ));
    if length(related_files) < 5
        include(file,1) = 0;
    else
        include(file,1) = 1;
    end
end
files = files(find(include),1);
% files = dir('*FTD*mSpikes_3.mat');  % where your .mat files are
% files = files(~contains({files.name}, 'ttx'));
% files = files(~contains({files.name}, 'TTX'));
% files = files(~contains({files.name}, 'DIV5'));
% files = files(~contains({files.name}, 'DIV4'));
% files = files(~contains({files.name}, 'DIV07'));
% files = files(~contains({files.name}, 'cleaned'));
% files = files(121:end);

sampling_fr = 25000;
sync_win_s = 0.01;

spikes_only         = 0;
elec_burst_stats    = 0;
Netw_burst_stats    = 0;
cor_ctrl            = 0; % 1 means correlate randomly shuffled trains; be sure to change the sync window for the ctrl
fc_figures          = 0;%change to 1 to add plots
g_figures           = 0;%graph theory
g_measures          = 1;

%% Mecp2 specific
% xlsname = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\XLS_unsorted\Mecp2 Culture Tracker 2022.xlsx';
% xlsname = ['C:\Users\alexd\OneDrive - University of Cambridge\',...
% 'Cam\PhD\Project\Data&CodeBackup\Data\XLS_unsorted\Mecp2 Culture Tracker 2022.xlsx'];
% sheet = 'MPA';
% [num,txt,raw] = xlsread(xlsname,sheet) ;
% genotypes_table = cell2table(raw(2:size(txt,1),1:size(txt,2)),'VariableNames',raw(1,1:size(txt,1)));
% clear num txt raw
% genotype_labels   = {'Wildtype','Mecp2^{-/Y}','Mecp2^{-/+}'};
% line_styles       = {'-',':'};
% DIVs        = {'DIV14','DIV21','DIV28','DIV35'};
% condition_labels  = {'Light on','Light off'};
% genotypes   = {'WT','KO','HE'};
% DIVs        = {14,21,28,35};
% conditions  = {0,1};
% gtcolours   = {[0 0 0],[1 0 0],[0 0 1]}

%% loop through
progressbar
for file = 1:length(files)
    %% load and get basic firing
    tic
    data=load(files(file).name,'*Spikes','channels','thresholds'); %may need to at try;catch;end loop if channels var not saved with spikes; could load from mat files
    try
        Spikes=data.mSpikes;
        disp('loaded mSpikes')
    catch
        Spikes=data.cSpikes;
        disp('loaded cSpikes')
    end
    
    try
        channels=data.channels;
    catch
        channels=load(strcat(files(file).name(1:18),'.mat'),'channels'); %some spike mats dont have channels var saved
        channels=channels.channels;
        disp('channels loaded from dat file')
    end
    
    try
        thresholds = data.thresholds;
        disp('loaded thresholds')
    catch
    end
    clear data
    % save space by making spike mat a logical array
    spikeMat=logical(full(Spikes));
    clear Spikes
    spikeCounts=sum(spikeMat);
    %remove ref channel spikes:
    spikeCounts(find(channels == 15)) = 0;        %cell containing n. spikes for each channel
    active_chanIndex=find(spikeCounts / (length(spikeMat)/sampling_fr) >= 0.01);
    ActiveSpikeCounts=spikeCounts(active_chanIndex);  %spikes of only active channels ('active'= >9)
    
    output(file).rec = files(file).name(1:end-4); %file name
    
    % find group info for organoid data
    %{ 
    if ~isempty(strfind(files(file).name(1:end-4),'Grp'))
        output(file).grp = output(file).rec(strfind(files(file).name(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(files(file).name(1:end-4),'Group'))
        output(file).grp = output(file).rec(strfind(files(file).name(1:end-4),'Group')+5);
    else
        output(file).grp = ['E'];
    end
    %}
    
    % for nanobody recs
    %{
    control_conditions = {  'BASELINE','BUFFER1_PRE','BUFFER1_POST','BUFFER2_PRE','BUFFER2_POST',...
        'BUFFER3_PRE','BUFFER3_POST','BUFFER4_PRE','BUFFER4_POST','FOLLOWUP'};
    nanobod_conditions = {  'BASELINE','PRE_0P3','POST_0P3','PRE_1P0','POST_1P0',...
        'PRE_3P0','POST_3P0','PRE_10','POST_10','FOLLOWUP'};
    condition_str_idx  = [strfind(files(file).name,'_DIV')+7 : strfind(files(file).name,'Spikes')-3 ];
    
    output(file).condition = files(file).name(condition_str_idx);
    
    related_files = files(~~contains({files.name}, files(file).name(1 : strfind(files(file).name,'_DIV'))  ));
    
    if ~isempty(related_files(~~contains({related_files.name}, 'BUFFER')))
        output(file).grp = ['Control'];
        output(file).grp_num = 1;
        output(file).condition_num = find(contains(control_conditions,output(file).condition));
    elseif  ~isempty(related_files(~~contains({related_files.name}, 'NB130')))
        output(file).grp = ['NB130'];
        output(file).grp_num = 2;
        output(file).condition = output(file).condition(7:end);
        output(file).condition_num = find(contains(nanobod_conditions,output(file).condition));
    elseif  ~isempty(related_files(~~contains({related_files.name}, 'NB138')))
        output(file).grp = ['NB138'];
        output(file).grp_num = 3;
        output(file).condition = output(file).condition(7:end);
        output(file).condition_num = find(contains(nanobod_conditions,output(file).condition));
    end
    clear related_files
    %}
    
    % get mecp2 group and div
    % read in table to specific sheet for this culture type using first
    % three characters experiment code
    gtdata = readtable(['C:\Users\alexd\OneDrive - University of Cambridge\',...
        'Cam\PhD\Project\Data&CodeBackup\Data\XLS_unsorted\Mecp2 Culture Tracker 2022.xlsx'],...
        'sheet',files(file).name(1:3));    
    testfilename = files(file).name;
        if strcmp(testfilename(4),'_')
            testfilename(4)='';
        end
        expID = testfilename(1:9);
        pupID = str2num(testfilename(11 : strfind(files(file).name,'_DIV')-3));
        gtidx = find(strcmp(gtdata.CultureID,expID  ) &...
            gtdata.PupID==pupID            );
        div  = str2num(testfilename(strfind(files(file).name,'_DIV')+3:strfind(files(file).name,'_DIV')+4));
        
    output(file).Culture    = testfilename(1:strfind(files(file).name,'_DIV')-2);
    output(file).Age        = testfilename(strfind(files(file).name,'_DIV'):strfind(files(file).name,'_DIV')+4);
    output(file).DIV        = round( div / 7 ) * 7;
    output(file).Genotype   = gtdata.Genotype{gtidx};
    
    if   isempty(strfind(output(file).rec,'ttx')) & isempty(strfind(output(file).rec,'TTX'))
        
        output(file).ttx = 0; %1 means ttx %0 means no ttx
        
    elseif   ~isempty(strfind(output(file).rec,'ttx')) | ~isempty(strfind(output(file).rec,'TTX'))
        
        output(file).ttx = 1;
        
    else
        disp('error finding whether file involved ttx')
    end
    
    
    output(file).spikes = spikeCounts;        %cell containing n. spikes for each channel
    if ~~exist('thresholds')
        output(file).thresholds = thresholds;
    end
    
    FiringRates=ActiveSpikeCounts/(length(spikeMat)/sampling_fr); %firing rate in seconds
    
    %if length(ActiveSpikeCounts)<4
    %output(file).mean_FR = 0; %if there are fewer than 4 active channels, set to 0 and exclude culture
    %output(file).median_FR = 0;
    %else
    output(file).total_FR = sum(FiringRates);
    output(file).mean_FR = round(mean(FiringRates),3);
    output(file).CV_FR = round(std(FiringRates),3)/round(mean(FiringRates),3);
    output(file).skew_FR = round(skewness(FiringRates),3);
    output(file).sd_FR = round(std(FiringRates),3);
    output(file).sem_FR = round(std(FiringRates)/(sqrt(length(ActiveSpikeCounts))),3);
    output(file).median_FR = round(median(FiringRates),3);
    output(file).iqr_FR = round(iqr(FiringRates),3);
    output(file).mad_FR = mad(FiringRates,1); % median absolute deviation
    %end
    output(file).N_active_Es = length(ActiveSpikeCounts);
    %get rid of NaNs where there are no spikes; change to 0
    if isnan(output(file).mean_FR);
        output(file).mean_FR=0;
    else
    end
    if isnan(output(file).median_FR);
        output(file).median_FR=0;
    else
    end
    disp('spike stats done')
    
    %% get within elec burst info
    if spikes_only ~= 1
        if elec_burst_stats == 1
            N = 6;
            [burstData] = fcn_bursts_within_bakkum(spikeMat,N,sampling_fr);
            
            output(file).array_mean_ISInTh            = burstData.burst_detection.array_mean_ISInTh;
            output(file).array_std_ISInTh             = burstData.burst_detection.array_std_ISInTh;
            output(file).n_bursting_elecs             = burstData.burst_detection.n_bursting_elecs;
            output(file).array_mean_burstRates        = burstData.burst_summaries.array_mean_burstRates;
            output(file).array_std_burstRates         = burstData.burst_summaries.array_std_burstRates;
            output(file).array_mean_mean_burstLengths = burstData.burst_summaries.array_mean_mean_burstLengths;
            output(file).array_std_mean_burstLengths  = burstData.burst_summaries.array_std_mean_burstLengths;
            output(file).array_mean_std_burstLengths  = burstData.burst_summaries.array_mean_std_burstLengths;
            output(file).array_std_std_burstLengths   = burstData.burst_summaries.array_std_std_burstLengths;
            output(file).array_mean_mean_IBI          = burstData.burst_summaries.array_mean_mean_IBI;
            output(file).array_std_mean_IBI           = burstData.burst_summaries.array_std_mean_IBI;
            output(file).array_mean_std_IBI           = burstData.burst_summaries.array_mean_std_IBI;
            output(file).array_std_std_IBI            = burstData.burst_summaries.array_std_std_IBI;
            output(file).array_mean_cv_IBI            = burstData.burst_summaries.array_mean_cv_IBI;
            output(file).array_std_cv_IBI             = burstData.burst_summaries.array_std_cv_IBI;
            output(file).array_mean_fracInBursts      = burstData.burst_summaries.array_mean_fracInBursts;
            output(file).array_std_fracInBursts       = burstData.burst_summaries.array_std_fracInBursts;
            output(file).array_sem_fracInBursts       = burstData.burst_summaries.array_sem_fracInBursts;
            output(file).array_sem_cv_IBI             = burstData.burst_summaries.array_sem_cv_IBI;
            output(file).array_sem_std_IBI            = burstData.burst_summaries.array_sem_std_IBI;
            output(file).array_sem_std_burstLengths   = burstData.burst_summaries.array_sem_std_burstLengths;
            output(file).array_sem_mean_IBI           = burstData.burst_summaries.array_sem_mean_IBI;
            output(file).array_sem_mean_burstLengths  = burstData.burst_summaries.array_sem_mean_burstLengths;
            output(file).array_sem_burstRates         = burstData.burst_summaries.array_sem_burstRates;
            
            output(file).array_mean_within_FR         = burstData.within_burst_summaries.array_mean_mean_WithinBurstFR;
            output(file).array_mean_cv_ISI_w          = burstData.within_burst_summaries.array_mean_cv_ISI_w;
            output(file).array_mean_mean_ISI_w        = burstData.within_burst_summaries.array_mean_mean_ISI_w;
            
            output(file).array_mean_outside_FR        = burstData.outside_burst_summaries.array_mean_mean_OutsideBurstFR;
            output(file).array_mean_cv_ISI_o          = burstData.outside_burst_summaries.array_mean_cv_ISI_o;
            output(file).array_mean_mean_ISI_o        = burstData.outside_burst_summaries.array_mean_mean_ISI_o;
        end
    end
    
    %% get network burst info
    if spikes_only ~= 1
        if Netw_burst_stats == 1
            %turn warning off
            warning('off','MATLAB:nearlySingularMatrix');
            
            samplingRate=sampling_fr;
            method ='Bakkum';
            %note, Set N = 30 (min number of bursts)
            %ensure bursts are excluded if fewer than 3 channels (see inside burst Detect
            %function)
            %to change min channels change line 207 of burstDetect.m
            %to change N (min n spikes) see line 170 of burstDetect.m
            % N = 30; minChan = 5;
            N = 40; minChan = 3;
            %{
            try
                [burstMatrix, burstTimes, burstChannels, ISInTh] = burstDetect(spikeMat, method, samplingRate,N,minChan);
                nBursts=length(burstMatrix);
                %trainCombine=sum(spikeMatrix, 2);
                if length(burstMatrix)>0
                    for Bst=1:length(burstMatrix)
                        sp_in_bst(Bst)=sum(sum(burstMatrix{Bst,1}));
                        train = sum(burstMatrix{Bst,1},2);%sum across channels
                        train(find(train>1))=1; %re-binarise
                        sp_times = find(train==1);
                        sp_times2= sp_times(2:end);
                        ISI_within = sp_times2 - sp_times(1:end-1);
                        mean_ISI_w(Bst) = round(mean(ISI_within)/sampling_fr,3); %in s with 3 d.p.
                        std_ISI_w(Bst) = round(std(ISI_within)/sampling_fr,3);
                        chans_involved(Bst) = length(burstChannels{Bst,1});
                        
                        NBLength(Bst) = size(burstMatrix{Bst,1},1)/sampling_fr;
                        
                        clear ISI_within sp_times sp_times2 train
                    end
                    % sp_in_bst=sum(sp_in_bst);
                else
                    disp('no bursts detected')
                    sp_in_bst=0;
                end
                
                train = sum(spikeMat,2);%sum across channels
                train(find(train>1))=1; %re-binarise
                sp_times = find(train==1);
                sp_times2= sp_times(2:end);
                ISI_all = sp_times2 - sp_times(1:end-1);
                
                clear train
                
                %CV of IBI
                %CV = st dev / mean
                %get IBIs
                end_times = burstTimes(1:end-1,2); %-1 because no IBI after end of last burst
                sta_times = burstTimes(2:end,1); %start from burst start time 2
                IBIs      = sta_times -end_times;
                % calculate CV of IBI and non need to convert from samples to seconds
                % (as relative measure it would be the same)
                
                % get ISIs outside bursts
                if burstTimes ~= 0
                    outside_burstStartTimes  = burstTimes(1:end-1,2);
                    outside_burstEndTimes    = burstTimes(2:end,1);
                    outside_burstTimes       = [outside_burstStartTimes, outside_burstEndTimes];
                    
                    % Add period before first burst unless the recording starts with a burst
                    if burstTimes(1,1) > 1
                        outside_burstTimes = [1 burstTimes(1,1) ; outside_burstTimes];
                    end
                    
                    % Add period after the last burst, unless the recording ends with a burst
                    if burstTimes(end,2) < length(spikeMat)
                        outside_burstTimes = [outside_burstTimes ; burstTimes(end,2) length(spikeMat)];
                    end
                    
                    % need to add option to save outside netw burst and
                    % inside netw burst matrix
                    outside_burstMatrix = sparse(zeros(size(spikeMat)));
                    for nb = 1:length(outside_burstTimes)
                        nb;
                        outside_burstMatrix(outside_burstTimes(nb,1)+1:outside_burstTimes(nb,2)-1,:) = spikeMat(outside_burstTimes(nb,1)+1:outside_burstTimes(nb,2)-1,:);
                    end
                    inside_burstMatrix  = sparse(zeros(size(spikeMat)));
                    for nb = 1:length(burstTimes)
                        inside_burstMatrix(burstTimes(nb,1):burstTimes(nb,2),:) = spikeMat(burstTimes(nb,1):burstTimes(nb,2),:);
                    end
                    
                    %                     % plot Netw burst
                    %                     MEAgraphics(1)
                    %                     tiledlayout(3,1)
                    %                     nexttile
                    %                     % imagesc(log10(downSampleSum(full(spikeMat),300))');            xticklabels([]);    ylabel('Electrode #');
                    %                     imagesc(downSampleSum(full(spikeMat),300)');            xticklabels([]);    ylabel('Electrode #');
                    %                     title('All spikes','FontWeight','Normal');
                    %                     ax = gca;
                    %                     ax.TitleHorizontalAlignment = 'left';
                    %                     colorbar('Visible','off');
                    %                     % caxis(log10([0.01 100]));
                    %                     caxis([0 100]);
                    %                     aesthetics
                    %
                    %                     nexttile
                    %                     % imagesc(log10(downSampleSum(full(outside_burstMatrix),300))'); xticklabels([]);    ylabel('Electrode #');
                    %                     imagesc(downSampleSum(full(outside_burstMatrix),300)'); xticklabels([]);    ylabel('Electrode #');
                    %                     title('Spikes outside network bursts','FontWeight','Normal');
                    %                     ax = gca;
                    %                     ax.TitleHorizontalAlignment = 'left';
                    %                     cb = colorbar; cb.Label.String = 'Spikes/s';
                    %                     % caxis(log10([0.01 100]));
                    %                     caxis([0 100]);
                    %                     aesthetics
                    %
                    %                     nexttile
                    %                     % imagesc(log10(downSampleSum(full(inside_burstMatrix),300))');  xlabel('Time (s)'); ylabel('Electrode #');
                    %                     imagesc(downSampleSum(full(inside_burstMatrix),300)');  xlabel('Time (s)'); ylabel('Electrode #');
                    %                     title('Spikes within network bursts','FontWeight','Normal');
                    %                     ax = gca;
                    %                     ax.TitleHorizontalAlignment = 'left';
                    %                     colorbar('Visible','off');
                    %                     % caxis(log10([0.01 100]));
                    %                     caxis([0 100]);
                    %                     aesthetics
                    
                else
                    disp('no network bursts detected')
                end
                
                % spikes in bursts, out of bursts and in total for all
                % elecs
                s_in_bursts     = full(sum(inside_burstMatrix));
                s_out_bursts    = full(sum(outside_burstMatrix));
                s_total         = sum(spikeMat);
                
                % get ISIs inside and outside bursts
                for elec = 1:size(spikeMat,2)
                    sp_times    = find(outside_burstMatrix(:,elec)==1);
                    sp_times2   = sp_times(2:end);
                    ISI_out     = sp_times2 - sp_times(1:end-1);
                    mean_ISI_out_elec(elec) = nanmean(ISI_out/sampling_fr);
                    std_ISI_out_elec(elec)  = nanstd(ISI_out/sampling_fr);
                    
                    clear sp_times sp_times2
                    
                    sp_times    = find(inside_burstMatrix(:,elec)==1);
                    sp_times2   = sp_times(2:end);
                    ISI_in      = sp_times2 - sp_times(1:end-1);
                    mean_ISI_in_elec(elec)  = nanmean(ISI_in/sampling_fr);
                    std_ISI_in_elec(elec)   = nanstd(ISI_in/sampling_fr);
                    
                    clear sp_times sp_times2 ISI_out ISI_in
                end
                
                % NOTE: these are based on the ISI across all channels!!!
                output(file).ISIn_NetwBurst_Threshold = ISInTh;
                output(file).num_Netw_bursts = length(burstMatrix);
                output(file).mean_count_chans_involved_in_Netw_bursts = mean(chans_involved);
                output(file).std_count_chans_involved_in_Netw_bursts = std(chans_involved);
                output(file).mean_Netw_Bst_length_s = mean(NBLength);
                output(file).std_Netw_Bst_length_s = std(NBLength);
                % note mean_num_spikes_inNetwBursts is the number of spikes
                % that occurs in the average burst i.e. how many spikes are
                % there in a burst... not the total number of spikes in all
                % bursts
                output(file).mean_num_spikes_inNetwBursts = mean(sp_in_bst);
                output(file).std_num_spikes_inNetwBursts = std(sp_in_bst);
                output(file).median_num_spikes_inNetwBursts = median(sp_in_bst);
                output(file).iqr_num_spikes_inNetwBursts = iqr(sp_in_bst);
                output(file).mean_proportion_chans_involved_in_Netw_bursts = mean(chans_involved)/length(active_chanIndex);
                output(file).std_proportion_chans_involved_in_Netw_bursts = std(chans_involved./length(active_chanIndex));
                
                % output(file).mean_mean_ISI_withinNetw_bursts_ms  = mean(mean_ISI_w);
                % output(file).std_mean_ISI_withinNetw_bursts_ms  = std(mean_ISI_w);
                % output(file).mean_std_ISI_withinNetw_bursts_ms  = mean(std_ISI_w);
                % output(file).std_std_ISI_withinNetw_bursts_ms  = std(std_ISI_w);
                output(file).mean_mean_ISI_withinNetw_bursts_ms     = nanmean(mean_ISI_in_elec);
                output(file).std_mean_ISI_withinNetw_bursts_ms      = nanstd(mean_ISI_in_elec);
                output(file).mean_std_ISI_withinNetw_bursts_ms      = nanmean(std_ISI_in_elec);
                output(file).std_std_ISI_withinNetw_bursts_ms       = nanstd(std_ISI_in_elec);
                
                output(file).mean_mean_ISI_outsideNetw_bursts_ms     = nanmean(mean_ISI_out_elec);
                output(file).std_mean_ISI_outsideNetw_bursts_ms      = nanstd(mean_ISI_out_elec);
                output(file).mean_std_ISI_outsideNetw_bursts_ms      = nanmean(std_ISI_out_elec);
                output(file).std_std_ISI_outsideNetw_bursts_ms       = nanstd(std_ISI_out_elec);
                
                CV_of_ISI_in_NetwBursts_elec    = std_ISI_in_elec./mean_ISI_in_elec;
                CV_of_ISI_out_NetwBursts_elec   = std_ISI_out_elec./mean_ISI_out_elec;
                
                output(file).mean_CV_of_ISI_in_NetwBursts      = nanmean(CV_of_ISI_in_NetwBursts_elec);
                output(file).std_CV_of_ISI_in_NetwBursts       = nanstd(CV_of_ISI_in_NetwBursts_elec);
                output(file).mean_CV_of_ISI_out_NetwBursts     = nanmean(CV_of_ISI_out_NetwBursts_elec);
                output(file).std_CV_of_ISI_out_NetwBursts      = nanstd(CV_of_ISI_out_NetwBursts_elec);
                
                % interspike intervals across entire collapsed spike train
                % output(file).mean_ISI_all = round(mean(ISI_all)/sampling_fr,3);
                % output(file).std_ISI_all = round(std(ISI_all)/sampling_fr,3);
                
                output(file).CVofINetwBI    = round((std(IBIs/sampling_fr)/mean(IBIs/sampling_fr)),3); %3 decimal places
                output(file).mean_IBI       = round(mean(IBIs/sampling_fr),3);
                output(file).std_IBI        = round(std(IBIs/sampling_fr),3);
                % clear unneeded vars
                clear end_times sta_times IBIs
                output(file).Netw_BurstRate=round(60*(nBursts/(length(spikeMat(:,1))/samplingRate)),3);
                
                % output(file).mean_frac_in_Netw_burst = round(sum(sp_in_bst)/sum(sum(spikeMat)),3);
                % output(file).std_frac_in_Netw_burst = round(sum(sp_in_bst)/sum(spikeMat),3);
                output(file).mean_frac_in_Netw_burst    = round(nanmean(s_in_bursts./s_total),3);
                output(file).std_frac_in_Netw_burst     = round(nanstd(s_in_bursts./s_total),3);
                % sanity check that proportion in + out = total ... works.
                % when making outburst matrix I ensured outbursttimes did
                % not overlap with bursttimes
                % [s_out_bursts'+s_in_bursts' s_total']
                % output(file).mean_frac_out_Netw_burst    = round(nanmean(s_out_bursts./s_total),3);
                % output(file).std_frac_out_Netw_burst     = round(nanstd(s_out_bursts./s_total),3);
                
                
                
                %need to go intro burst detect and edit as it is not deleting the bursts
                %with <5 channels from burstChannels and burstTimes hence they are longer
                %need this for easier plotting of burst
                disp('burst stats done')
            catch
                disp('burst stats failed')
            end
            %}
            
            % updated burst section
            try
                % may require downsampling
                % updated nburst function to keep inside and outside burst
                % matrices as sparse logical arrays so downsampling no
                % longer needed
                down_fr = samplingRate; % samples /s; 1000 hz is not enough resolution; need at least 5000
                downFactor = samplingRate/down_fr;
                downSpikeMatrix = sparse(downSampleSum(spikeMat,length(spikeMat)/downFactor));
                downSpikeMatrix(downSpikeMatrix>1)=1;

                [network_burstData] = fcn_bursts_network_bakkum(downSpikeMatrix,N,down_fr,minChan);
                clear downSpikeMatrix
                
                output(file).NetwBurst_ISIn_Threshold       = network_burstData.burst_detection.array_mean_ISInTh;
                output(file).NetwBurst_num_elecs_involved   = network_burstData.burst_detection.mean_n_bursting_elecs;
                output(file).NetwBurst_CV_elecs_involved    = network_burstData.burst_detection.std_n_bursting_elecs/network_burstData.burst_detection.mean_n_bursting_elecs;
                
                output(file).NetwBurst_per_min              = network_burstData.burst_summaries.array_mean_burstRates;
                % output(file).num_Netw_bursts = network_burstData.burst_summaries.array_mean_burstRates;
                output(file).NetwBurst_CV_IBI               = network_burstData.burst_summaries.array_mean_cv_IBI
                output(file).NetwBurst_Frac_in_bursts       = network_burstData.burst_summaries.array_mean_fracInBursts
                output(file).NetwBurst_Length               = network_burstData.burst_summaries.array_mean_mean_burstLengths
                
                output(file).NetwBurst_CV_ISI_within        = network_burstData.within_burst_summaries.array_mean_cv_ISI_w;
                output(file).NetwBurst_FR_within            = network_burstData.within_burst_summaries.array_mean_mean_WithinBurstFR;
                output(file).NetwBurst_ISI_within           = network_burstData.within_burst_summaries.array_mean_mean_ISI_w;
                
                output(file).NetwBurst_CV_ISI_outside       = network_burstData.outside_burst_summaries.array_mean_cv_ISI_o;
                output(file).NetwBurst_FR_outside           = network_burstData.outside_burst_summaries.array_mean_mean_OutsideBurstFR;
                output(file).NetwBurst_ISI_outside          = network_burstData.outside_burst_summaries.array_mean_mean_ISI_o;

                output(file).NetwBurst_CV_Bst_Lengths       = network_burstData.burst_summaries.array_mean_std_burstLengths/network_burstData.burst_summaries.array_mean_mean_burstLengths;
                
            catch
                disp('burst stats failed')
            end
        end
        %% functional connectivity
        method = 'tileCoef';
        %get adjM - ideally already created by batch get adjM to save time when re
        %analysing
        
        if length(active_chanIndex)<3 %set all output to 0 etc. if no activity
            %else run through analyses
            
            output(file).meanSTTC=0;
            output(file).STTC_RCtrl=0;
            %output(file).netw_density=0;
            %output(file).meanDegree =0;
            %output(file).CC=0;
            %output(file).PL=58;
            %output(file).SW=0;
            %output(file).RC=0;
            
            disp('no active channels')
            
            
        else
            
            active_spike_mat = spikeMat(:,active_chanIndex);
            sync_win_s
            if ~exist(strcat(files(file).name(1:end-4), '_adjM_',num2str(sync_win_s), '.mat'))
                disp('could not find any adjM saved; update section to downsample and get adj matrix')
                % below four lines need updating to quickly get adjM if not
                % exist already
                active_adjM=getAdjM(active_spike_mat, method, 0,sync_win_s);
                adjM=active_adjM;
                active_adjM=active_adjM-eye(size(active_adjM));
                active_adjM = weight_conversion(active_adjM,'autofix');
                % active_adjM(find(isnan(active_adjM)))=0;
                
                %adjM = getAdjM(spikeMat, method, 0,0.05); %0 means no downsampling; 0.05 is sync window in s
                
                
            else
                disp('loading adjM')
                adjM=load(strcat(files(file).name(1:end-4), '_adjM_',num2str(sync_win_s), '.mat'));
                adjM=adjM.adjM;
                adjM1= adjM-eye(size(adjM));
                adjM1(find(isnan(adjM1)))=0;
                active_adjM=adjM1(active_chanIndex,active_chanIndex); %for extracting
                active_adjM = weight_conversion(active_adjM,'autofix');
                %part of full adjmat but now just compute active adjmat
                
                
            end
            
            %mean correlation of active channels
            % output(file).meanSTTC= round(sum(sum(active_adjM))/(length(active_adjM)*(length(active_adjM)-1)),3);
            % output(file).semSTTC = round(std(nonzeros(triu(active_adjM)))/(sqrt(length(nonzeros(triu(active_adjM))))),3);
            % output(file).sdSTTC = round(std(nonzeros(triu(active_adjM))),3);
            strengths       = strengths_und(active_adjM)./length(active_adjM);
            sttc_sd_elecs   = std(active_adjM);
            
            output(file).mean_mean_STTC     = round(mean(strengths),3);
            output(file).sem_mean_STTC      = round(std(strengths)/(sqrt(length(strengths))),3);
            output(file).sd_mean_STTC       = round(std(strengths),3);
            output(file).cv_mean_STTC       = round(std(strengths),3) / round(mean(strengths),3);
            
            output(file).mean_sd_STTC     = round(mean(sttc_sd_elecs),3);
            output(file).sem_sd_STTC      = round(std(sttc_sd_elecs)/(sqrt(length(sttc_sd_elecs))),3);
            output(file).sd_sd_STTC       = round(std(sttc_sd_elecs),3);
            output(file).mean_cv_STTC     = round(mean(sttc_sd_elecs./strengths),3);
            
            weights = active_adjM;
            weights = weights - eye(size(weights));
            weights(find(isnan(weights))) = 0;
            weights(find(weights < 0)) = 0;
            mweights = mean(weights)';
            output(file).skewSTTC = skewness(mweights);
            
            %% fc control - randomise FR of each active electrode and compute STTC
            %create randomised controls for each recording and save their control
            %adjMs using the batch get adjM
            
            %load spike mat and find active channels
            
            %save memory; clear some vars
            %clear spikeMat
            clear cSpikes mSpikes
            clear burstMatrix burstChannels burstTimes Bst nBursts sp_in_bst
            
            if cor_ctrl == 1
                
                
                %                 bin_size = 0.01;
                %                 length(spikeMat(:,1))/25;
                num_iterations = 100;
                for k = 1:num_iterations %for running multiple randomisation iterations
                    disp(strcat(num2str(1+num_iterations-k),'_iterations_remaining'))
                    tic
                    for elec = 1:length(spikeMat(1,:))
                        %{
 approach:
                take each 1 ms window and randomly shuffle them
                get a vector of indices for each 1 ms period
                randomly shuffle this
                then resort the spike train using this index
                        %}
                        
                        down_factor = 25;
                        fs = 25000;
                        fr = sum(spikeMat(:,elec))/(length(spikeMat(:,elec))/fs); %firing rate %divide by seconds
                        tSim = length(spikeMat(:,elec))/fs; %duration of simulation in s
                        nTrials = 1;%number of trials
                        dt = (1/fs)*down_factor; %time bins; set to 1 ms currently
                        [rand_spike_vec, tVec] = poissonSpikeGen(fr, tSim, nTrials,dt);
                        rand_spikeMat(elec,:) = abs(rand_spike_vec);
                    end
                    rand_spikeMat = rand_spikeMat';
                    
                    %plot to compare spike mat spike counts and rand spike
                    %counts over x iterations
                    % a(k,:) = sum(spikeMat) - sum(rand_spikeMat);
                    
                    %get ctrl adjM
                    
                    rec_length_samps = length(rand_spikeMat);
                    rec_length_s = rec_length_samps/(sampling_fr/down_factor);
                    num_samples = 1000 * rec_length_s; %defaults to 1 sample per ms
                    ds_rate = length(rand_spikeMat)/length(spikeMat); %scalar for time_window %downsampling rate
                    sync_win = sync_win_s * ds_rate;
                    adjM2 = getAdjM(rand_spikeMat, method, 0,sync_win); %0 means no downsampling; 0.05 is sync window in s
                    adjM2 = adjM2(active_chanIndex,active_chanIndex);
                    adjM2 = adjM2 - eye(size(adjM2));
                    adjM2(find(isnan(adjM2)))=0;
                    %save every random adjM
                    adjM2s(:,:,k) = adjM2;
                    %save each mean STTC for the randoms
                    randSTTC(k) = round(sum(sum(adjM2))/(length(adjM2)*(length(adjM2)-1)),3);
                    skewSTTC(k) = skewness(sum(adjM2) / (length(adjM2)-1));
                    %take one of the random spike mats (not enough memory
                    %to take all)
                    if k == 1
                        rand_spikeMat_sparse = sparse(rand_spikeMat);
                    end
                    
                    clear adjM2 rand_spikeMat
                    toc
                    
                end %for running multiple randomisation iterations
                %                    adjM2 = mean(adjM2s,3);%get mean of random adjMs
                adjM2 = adjM2s(:,:,1);%get first random adjM
                
                % save first control adjM
                %                 if file == 1
                %                     fileNameSpikes = strcat(files(file).name(1:end-4), '_CTRL_spikeMat.mat');
                %                     fileNameadjM = strcat(files(file).name(1:end-4), '_CTRL_adjM_3.mat');
                %                     save(fileNameSpikes, 'rand_spikeMat_sparse','-v7.3');
                %                     save(fileNameadjM, 'adjM2','-v7.3');
                %                     disp('saved first ctrl mat')
                %                 end
                
                clear adjM2 rand_spikeMat
                
                output(file).STTC_RCtrl_mean = mean(randSTTC);
                output(file).STTC_RCtrl_skew = mean(skewSTTC);
                output(file).nrmlsd_cnnctvty = output(file).meanSTTC/output(file).STTC_RCtrl_mean;
                %figure; plot(mean(a))
                
                if fc_figures == 1
                    %compare
                    figure;
                    imagesc(active_adjM,[min(min(active_adjM)),1]); %set colouring limits!
                    xlabel('Channel ID');
                    ylabel('Channel ID');
                    aesthetics
                    box off
                    set(gca,'TickLength',[0 0])
                    %set(gca,'xtick',[])
                    xticks(linspace(1,length(active_adjM),length(active_adjM))) %adjust to n active channels
                    yticks(linspace(1,length(active_adjM),length(active_adjM))) %adjust to n active channels
                    xticklabels(channels(active_chanIndex))
                    yticklabels(channels(active_chanIndex))
                    cbFC=colorbar
                    cbFC.Limits=[min(min(active_adjM)),1]
                    
                    
                    figure;
                    imagesc(adjM2,[min(min(active_adjM)),1])
                    aesthetics
                    box off
                    set(gca,'TickLength',[0 0])
                    xticks(linspace(1,length(active_adjM),length(active_adjM))) %adjust to n active channels
                    yticks(linspace(1,length(active_adjM),length(active_adjM))) %adjust to n active channels
                    xticklabels(channels(active_chanIndex))
                    yticklabels(channels(active_chanIndex))
                    cbFD=colorbar
                    cbFD.Limits=[min(min(active_adjM)),1]
                else
                end
            end
            disp('fc stats done')
            %% Graph metrics
            if g_measures == 1
                disp('beginning graph measures')
                % Basic: calculate over three absolute thresholds
                % 0.3 0.5 0.7
                count1 = 1; %to track threshold iterations
                for cutoff = [linspace(5, 95, 19) / 100]
                    
                    threshold = cutoff;
                    
                    % % choose either whole mat or only active channels
                    %   edges=adjM1;
                    edges=active_adjM;
                    
                    edges(edges < threshold) = 0.0001;
                    
                    buEdges=edges;
                    buEdges(find(buEdges==0.0001))=0;
                    buEdges(find(buEdges>0.0001))=1;
                    %remove nodes with no edges in binary network
                    active_nodes_vec=find(sum(buEdges)>0);
                    buEdges=buEdges(active_nodes_vec,active_nodes_vec);
                    ActiveChanLabels=channels(active_chanIndex);
                    NodeLabels=ActiveChanLabels(active_nodes_vec);
                    
                    if length(buEdges)==0 %if there are no nodes with edges above threshold
                        netw_density(count1) = 0;%else it would be NaN with below formula
                    else
                        netw_density(count1) = sum(sum(buEdges))/(length(buEdges)*(length(buEdges)-1));
                    end
                    meanDegree(count1)       = mean(sum(buEdges));
                    netw_size (count1)       = length(buEdges); %num. nodes
                    count1 = count1 + 1;
                end
                
                output(file).netw_density   =   mean(netw_density);
                output(file).meanDegree     =   mean(meanDegree);
                output(file).skewDegree     =   skewness(meanDegree);
                output(file).netw_size      =   mean(netw_size);
                
                disp('simple graph measures done, starting complex measures')
                % Complex: calculate over three proportion thresholds
                % 40, 60 and 80 %
                
                RC_score = zeros(length(channels),1);
                
                count2 = 1;
                for cutoff = [60]
                    
                    % get to 40% of connections
                    % threshold = prctile(adjM(:), cutoff);
                    adjM_th_1 = adjM - 6 * eye(size(adjM));
                    w = find(adjM_th_1 >-5);
                    threshold = prctile(adjM_th_1(w),cutoff);
                    clear adjM_th_1 w
                    
                    % %                     choose either whole mat or only active channels
                    %                     edges=adjM1;
                    edges=active_adjM;
                    
                    %                     edges(edges < threshold) = 0.0001;
                    
                    edges(edges >= threshold) = 1;
                    edges(edges < threshold)  = 0;
                    edges = weight_conversion(edges, 'autofix');
                    buEdges=edges;
                    %                     buEdges(find(buEdges==0.0001))=0;
                    %                     buEdges(find(buEdges>0.0001))=1;
                    % correction for remaining values
                    %                     buEdges(find(buEdges<0.0001))=0;
                    %remove nodes with no edges in binary network
                    active_nodes_vec=find(sum(buEdges)>0);
                    buEdges=buEdges(active_nodes_vec,active_nodes_vec);
                    ActiveChanLabels=channels(active_chanIndex);
                    NodeLabels=ActiveChanLabels(active_nodes_vec);
                    
                    
                    
                    if round(max(sum(buEdges))/5)*5 > 0 & length(buEdges) >= 5 % check for enough nodes
                        % clustering
                        C = clustering_coef_bu(buEdges);
                        C = mean(C);
                        % path length and global efficiency GE
                        D           = distance_bin(buEdges);
                        [PL, GE]    = charpath(D,0,0); % 0 and 0 excludes self connections and infinities (non connected) respectively
                        % small world
                        SW = C/PL;
                        % rich club and betweenness centrality
                        k               =   max(sum(buEdges));
                        [Rc,Nk,Ek]      =   rich_club_bu(buEdges,k);
                        RC              =   max(Rc);
                        maxKrand        =   min(find(Rc==RC));
                        if ~isempty(maxKrand)
                            RC_nodes_vec    =   find(sum(buEdges) >= maxKrand);
                            num_RC_nodes(count2)        =   length(RC_nodes_vec);
                            RC_node_IDs{count2,1}       =   NodeLabels(RC_nodes_vec);
                            current_RC_nodes = RC_node_IDs{count2,1};
                            for chan = 1:length(RC_node_IDs{count2,1})
                                RC_score(find(channels == current_RC_nodes(chan))) = RC_score(find(channels == current_RC_nodes(chan))) +1;
                            end
                            clear current_RC_nodes
                        end
                        
                        if ~isempty(maxKrand)
                            BC_vec          =   betweenness_bin(buEdges);
                            BC              =   mean(BC_vec(RC_nodes_vec));
                            BC_norm         =   BC / ((length(buEdges)-1) * (length(buEdges)-2));
                        end
                        
                        % normalisation
                        disp('computing normalisations')
                        p_iterations = 20;
                        for r_i = 1:p_iterations
                            
                            [R,eff]     = randmio_und(buEdges, 100); %preserves density and distribution
                            
                            Nnodes      = length(buEdges);
                            Nedges      = sum(sum(buEdges))/2; % divide by two for non directional network
                            R2          = makerandCIJ_und(Nnodes,Nedges);
                            
                            CrVEC(:,r_i)= clustering_coef_bu(R);
                            Cr(r_i)     = mean(CrVEC(:,r_i));
                            CrVEC2(:,r_i)= clustering_coef_bu(R2);
                            Cr2(r_i)     = mean(CrVEC2(:,r_i));
                            
                            Dr                      = distance_bin(R);
                            [PLr(r_i),GEr(r_i)]     = charpath(Dr,0,0); % 0 and 0 excludes self connections and infinities (non connected) respectively
                            Dr2                     = distance_bin(R2);
                            [PLr2(r_i),GEr2(r_i)]   = charpath(Dr2,0,0); % 0 and 0 excludes self connections and infinities (non connected) respectively
                            
                            SWr(r_i)         =  Cr(r_i) / PLr(r_i);
                            SWr2(r_i)        =  Cr2(r_i) / PLr2(r_i);
                            
                            % for RC p value, I bootstrap. I take the mean
                            % of ten iterations 100 times. Avoids getting
                            % occasional RC values of 1 when they're
                            % generally low around say .3. Then round to 1
                            % decimal place so .9 will round up if that's
                            % the average
                            if ~isempty(maxKrand)
                                for RC_i = 1:10
                                    [R_rc1,eff]          = randmio_und(buEdges, 100);
                                    R_rc2               = makerandCIJ_und(Nnodes,Nedges);
                                    
                                    [Rcr,Nk1,Ek1]       = rich_club_bu(R_rc1,maxKrand);
                                    Rcr(find(isnan(Rcr))) = 0; % correction in case rich club coeff is NaN
                                    RCr_RC_i(RC_i)      = max(Rcr);
                                    
                                    [Rcr2,Nk2,Ek2]      = rich_club_bu(R_rc2,maxKrand);
                                    Rcr2(find(isnan(Rcr2))) = 0; % correction in case rich club coeff is NaN
                                    RCr2_RC_i(RC_i)     = max(Rcr2);
                                end
                                RCr(r_i)    = round(mean(RCr_RC_i),1);
                                RCr2(r_i)   = round(mean(RCr2_RC_i),1);
                                
                                % number of edges required to be in random network's rich club:
                                Rcr2K = find(Rcr2==max(Rcr2));
                                RcrK = find(Rcr==max(Rcr));
                                % find nodes in rand netw with >k edges
                                RCr_nodes_vec   =   find(sum(R) >= RcrK(1));
                                RCr2_nodes_vec  =   find(sum(R2) >= Rcr2K(1));
                                
                                BCr_vec          =   betweenness_bin(R);
                                BCr(r_i)         =   mean(BCr_vec(RCr_nodes_vec));
                                BCr_norm(r_i)    =   BCr(r_i)/((length(buEdges)-1) * (length(buEdges)-2));
                                
                                BCr2_vec          =   betweenness_bin(R2);
                                BCr2(r_i)         =   mean(BCr2_vec(RCr2_nodes_vec));
                                BCr2_norm(r_i)    =   BCr2(r_i)/((length(buEdges)-1) * (length(buEdges)-2));
                                clear BCr_vec BCr2_vec CrVEC CrVEC2
                            end
                        end
                        if ~isempty(maxKrand)
                            BCr2_norm = BCr2_norm(~isnan(BCr2_norm));
                        end
                        disp('normalisation computed')
                        
                        all_CC_raw(count2)  =   C;
                        all_CC(count2)      =   C/mean(Cr);
                        all_CC_p(count2)    =   length(find(Cr>=C))/p_iterations;
                        all_CC2(count2)     =   C/mean(Cr2);
                        all_CC2_p(count2)   =   length(find(Cr2>=C))/p_iterations;
                        
                        % note here I am asking is the path length
                        % significantly higher in the real than the random
                        % networks...
                        all_PL_raw(count2)  =   PL;
                        all_PL(count2)      =   PL/mean(PLr);
                        all_PL2(count2)     =   PL/mean(PLr2);
                        all_PL_p(count2)    =   length(find(PLr>=PL))/p_iterations;
                        all_PL2_p(count2)   =   length(find(PLr2>=PL))/p_iterations;
                        
                        % note here I am asking is the global efficiency
                        % significantly worse in the real than the random
                        % networks... if not, then we can say they are as
                        % efficient at global comms as a random network
                        all_GE_raw(count2)  =   GE;
                        all_GE(count2)      =   GE/mean(GEr);
                        all_GE_p(count2)    =   length(find(GEr<=GE))/p_iterations;
                        all_GE2(count2)     =   GE/mean(GEr2);
                        all_GE2_p(count2)   =   length(find(GEr2<=GE))/p_iterations;
                        
                        all_SW_raw(count2)  =   all_CC_raw(count2)/all_PL_raw(count2);
                        all_SW(count2)      =   all_CC(count2)/all_PL(count2);
                        all_SW2(count2)     =   all_CC2(count2)/all_PL2(count2);
                        all_SW_p(count2)    =   length(find(SWr>=SW))/p_iterations;
                        all_SW2_p(count2)   =   length(find(SWr2>=SW))/p_iterations;
                        
                        if ~isempty(maxKrand)
                            
                            all_RC_raw(count2)  =   RC;
                            all_RC(count2)      =   RC/mean(RCr);
                            all_RC2(count2)     =   RC/mean(RCr2);
                            all_RC_p(count2)    =   length(find(RCr>=RC))/p_iterations;
                            all_RC2_p(count2)   =   length(find(RCr2>=RC))/p_iterations;
                            
                            all_BC_raw(count2)  =   BC;
                            all_BC(count2)      =   BC/mean(BCr(~isnan(BCr)));
                            all_BC_p(count2)    =   length(find(BCr>=BC))/p_iterations;
                            all_BC2(count2)     =   BC/mean(BCr2(~isnan(BCr2)));
                            all_BC2_p(count2)   =   length(find(BCr2>=BC))/p_iterations;
                            
                            all_BC_n_raw(count2)=   BC_norm;
                            all_BC_n(count2)    =   BC_norm/mean(BCr_norm);
                            all_BC_n_p(count2)  =   length(find(BCr_norm>=BC_norm))/p_iterations;
                            all_BC2_n(count2)   =   BC_norm/mean(BCr2_norm);
                            all_BC2_n_p(count2) =   length(find(BCr2_norm>=BC_norm))/p_iterations;
                        else
                            
                            all_RC_raw(count2)  =    0;
                            all_RC(count2)      =    0;
                            all_RC2(count2)     =    0;
                            all_RC_p(count2)    =    1;
                            all_RC2_p(count2)   =    1;
                            
                            all_BC_raw(count2)  =    0;
                            all_BC(count2)      =    0;
                            all_BC2(count2)     =    0;
                            all_BC_p(count2)    =    1;
                            all_BC2_p(count2)   =    1;
                            
                            all_BC_n_raw(count2)=    0;
                            all_BC_n(count2)    =    0;
                            all_BC2_n(count2)   =    0;
                            all_BC_n_p(count2)  =    1;
                            all_BC2_n_p(count2) =    1;
                            
                        end
                    else
                        
                        all_CC_raw(count2)  =    0;
                        all_CC(count2)      =    0;
                        all_CC2(count2)     =    0;
                        all_CC_p(count2)    =    1;
                        all_CC2_p(count2)   =    1;
                        
                        all_PL_raw(count2)  =    length(active_chanIndex)-1;
                        all_PL(count2)      =    length(active_chanIndex)-1;
                        all_PL2(count2)     =    length(active_chanIndex)-1;
                        all_PL_p(count2)    =    1;
                        all_PL2_p(count2)   =    1;
                        
                        all_GE_raw(count2)  =    0;
                        all_GE(count2)      =    0;
                        all_GE2(count2)     =    0;
                        all_GE_p(count2)    =    1;
                        all_GE2_p(count2)   =    1;
                        
                        all_SW_raw(count2)  =    0;
                        all_SW(count2)      =    0;
                        all_SW2(count2)     =    0;
                        all_SW_p(count2)    =    1;
                        all_SW2_p(count2)   =    1;
                        
                        all_RC_raw(count2)  =    0;
                        all_RC(count2)      =    0;
                        all_RC2(count2)     =    0;
                        all_RC_p(count2)    =    1;
                        all_RC2_p(count2)   =    1;
                        
                        all_BC_raw(count2)  =    0;
                        all_BC(count2)      =    0;
                        all_BC2(count2)     =    0;
                        all_BC_p(count2)    =    1;
                        all_BC2_p(count2)   =    1;
                        
                        all_BC_n_raw(count2)=    0;
                        all_BC_n(count2)    =    0;
                        all_BC2_n(count2)   =    0;
                        all_BC_n_p(count2)  =    1;
                        all_BC2_n_p(count2) =    1;
                        
                    end
                    
                    count2 = count2 + 1;
                    disp('a prctle threshold iteration completed')
                    clear CrVEC CrVEC2
                    
                end
                RC_score = RC_score./(count2-1);
                
                if ~exist('num_RC_nodes')
                    num_RC_nodes = 0;
                    RC_node_IDs = {NaN};
                end
                
                output(file).num_RC_nodes_mean          = mean(num_RC_nodes);
                output(file).num_RC_nodes_maxthresh     = num_RC_nodes(count2-1);
                output(file).RC_node_IDs_allthresh      = channels(find(RC_score == 1));
                output(file).RC_node_IDs_maxthresh      = RC_node_IDs{count2-1,:};
                
                output(file).CCraw      =   mean(all_CC_raw);
                output(file).PLraw      =   mean(all_PL_raw);
                output(file).GEraw      =   mean(all_GE_raw);
                output(file).SWraw      =   mean(all_SW_raw);
                output(file).RCraw      =   mean(all_RC_raw);
                output(file).BCraw      =   mean(all_BC_raw);
                output(file).BC_nraw    =   mean(all_BC_n_raw);
                
                output(file).CCraw      =   all_CC_raw(count2-1);
                output(file).PLraw      =   all_PL_raw(count2-1);
                output(file).GEraw      =   all_GE_raw(count2-1);
                output(file).SWraw      =   all_SW_raw(count2-1);
                output(file).RCraw      =   all_RC_raw(count2-1);
                output(file).BCraw      =   all_BC_raw(count2-1);
                output(file).BC_nraw    =   all_BC_n_raw(count2-1);
                
                output(file).CC1        =   mean(all_CC);
                output(file).PL1        =   mean(all_PL);
                output(file).GE1        =   mean(all_GE);
                output(file).SW1        =   mean(all_SW);
                output(file).RC1        =   mean(all_RC);
                output(file).BC1        =   mean(all_BC);
                output(file).BC1_n      =   mean(all_BC_n);
                
                output(file).CC2        =   mean(all_CC2);
                output(file).PL2        =   mean(all_PL2);
                output(file).GE2        =   mean(all_GE2);
                output(file).SW2        =   mean(all_SW2);
                output(file).RC2        =   mean(all_RC2);
                output(file).BC2        =   mean(all_BC2);
                output(file).BC2_n      =   mean(all_BC2_n);
                
                output(file).CC1max_thresh     =   all_CC(count2-1);
                output(file).PL1max_thresh     =   all_PL(count2-1);
                output(file).GE1max_thresh     =   all_GE(count2-1);
                output(file).SW1max_thresh     =   all_SW(count2-1);
                output(file).RC1max_thresh     =   all_RC(count2-1);
                output(file).BC1max_thresh     =   all_BC(count2-1);
                output(file).BC1_nmax_thresh   =   all_BC_n(count2-1);
                
                output(file).CC2max_thresh     =   all_CC2(count2-1);
                output(file).PL2max_thresh     =   all_PL2(count2-1);
                output(file).GE2max_thresh     =   all_GE2(count2-1);
                output(file).SW2max_thresh     =   all_SW2(count2-1);
                output(file).RC2max_thresh     =   all_RC2(count2-1);
                output(file).BC2max_thresh     =   all_BC2(count2-1);
                output(file).BC2_nmax_thresh   =   all_BC2_n(count2-1);
                
                output(file).CC1_p        =   mean(all_CC_p);
                output(file).PL1_p        =   mean(all_PL_p);
                output(file).GE1_p        =   mean(all_GE_p);
                output(file).SW1_p        =   mean(all_SW_p);
                output(file).RC1_p        =   mean(all_RC_p);
                output(file).BC1_p        =   mean(all_BC_p);
                output(file).BC1_n_p      =   mean(all_BC_n_p);
                
                output(file).CC2_p        =   mean(all_CC2_p);
                output(file).PL2_p        =   mean(all_PL2_p);
                output(file).GE2_p        =   mean(all_GE2_p);
                output(file).SW2_p        =   mean(all_SW2_p);
                output(file).RC2_p        =   mean(all_RC2_p);
                output(file).BC2_p        =   mean(all_BC2_p);
                output(file).BC2_n_p      =   mean(all_BC2_n_p);
                
                output(file).CC1_pmax_thresh        =   all_CC_p(count2-1);
                output(file).PL1_pmax_thresh        =   all_PL_p(count2-1);
                output(file).GE1_pmax_thresh        =   all_GE_p(count2-1);
                output(file).SW1_pmax_thresh        =   all_SW_p(count2-1);
                output(file).RC1_pmax_thresh        =   all_RC_p(count2-1);
                output(file).BC1_pmax_thresh        =   all_BC_p(count2-1);
                output(file).BC1_n_pmax_thresh      =   all_BC_n_p(count2-1);
                
                output(file).CC2_pmax_thresh        =   all_CC2_p(count2-1);
                output(file).PL2_pmax_thresh        =   all_PL2_p(count2-1);
                output(file).GE2_pmax_thresh        =   all_GE2_p(count2-1);
                output(file).SW2_pmax_thresh        =   all_SW2_p(count2-1);
                output(file).RC2_pmax_thresh        =   all_RC2_p(count2-1);
                output(file).BC2_pmax_thresh        =   all_BC2_p(count2-1);
                output(file).BC2_n_pmax_thresh      =   all_BC2_n_p(count2-1);
                
                disp('complex graph measures complete')
                
            end
            
        end %end of loop if no active channels
    end
    disp('file done')
    toc
    progressbar(file/length(files));
    
    clearvars -except files file spikes_only elec_burst_stats Netw_burst_stats cor_ctrl ...
        fc_figures g_figures g_measures sampling_fr output sync_win sync_win_s gtdata
    
end


%% save
disp('saving...')
method = 'cSpikes';
threshold = '5SD';
fileName = strcat(method,'_' ,threshold, '_Mecp2_analysis_2022_Aug_10_ms_dt_STTC_graph_only', '.mat');
save(fileName, 'output');
xldata = struct2table(output);
xldata(:,'spikes') = [];
xldata(:,'thresholds') = [];
try
    xldata(:,'RC_node_IDs_allthresh') = [];
    xldata(:,'RC_node_IDs_maxthresh') = [];
catch
end
writetable(xldata, strcat(fileName(1:end-4),'.xlsx'))
