% check spikes for given set of recs ; 
% input = .mat data files

%cd to data files folder
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.1.MPT190515_mat_files'

%get spike matrices (if already done comment out)
%batchGetSpike

%if got spikes already go to spikes folder:
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\mats'

files = dir('*Spikes*.mat');  % where your .mat files are

sampling_fr = 25000;

for file = 1:length(files)
load(files(file).name,'*Spikes','channels');
Spikes=cSpikes;

spikeMat=full(Spikes);
spikeCounts=sum(spikeMat);
ActiveSpikeCounts=spikeCounts(find(spikeCounts>9));  %spikes of only active channels ('active'= >9)

output(file).rec = files(file).name(1:18); %file name
out_spikes=channels;
out_spikes(:,2)=sum(spikeMat)';
output(file).spikes = out_spikes;        %cell containing n. spikes for each channel

FiringRates=spikeCounts/(length(spikeMat)/sampling_fr); %firing rate in seconds
Active_FRs=FiringRates(find(spikeCounts>9)); %firing rates in Hz of active channels only

output(file).mean_FR = mean(Active_FRs);
output(file).median_FR = median(Active_FRs);
output(file).N_active_Es = length(ActiveSpikeCounts);

%get rid of NaNs where there are no spikes; change to 0
if isnan(output(file).mean_FR);
    output(file).mean_FR=0;
end
if isnan(output(file).median_FR);
    output(file).median_FR=0;
end


clear Spikes
clear spikeMat
clear spikeCounts
clear ActiveSpikeCounts
clear FiringRates
clear Active_FRs

progressbar(file/length(files));

end


        %% save 
        method = 'manuel';
        threshold = '7.5';
        fileName = strcat(method, threshold, '_stats', '.mat'); 
        save(fileName, 'output');
        