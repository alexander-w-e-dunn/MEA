%{
% Data structure
dm_spk_data = 
spike_times: {1×334 cell} % spike times in seconds (sampling rate is 20 kHz)
coordinates: [334×2 double] % in microns

Note changes since last analysis
- there are now no cells to remove (already removed this time)
- variable names for struct array , spike times and coordinates have
changed
- spike times are now in seconds not samples, so I have adjusted
downsampling
- negative relationships included as thresholding is run on abs(adjM)

as before
- coordinates are in microns
%}
clear all; close all

dataDir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\GNM_Data'; 
scriptsDir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
addpath(genpath(scriptsDir));
cd(dataDir);
files = dir('14*R.mat');

Hz=20000;
downHz = 1000; % one sample per ms
% downFactor = downHz/Hz; % needed for STTC
time_window_minutes = 10;
sample_window = time_window_minutes * 60 * downHz ; % in milliseconds

progressbar('files','cells')
for file = 1 : length(files)
    clearvars -except files downHz downFactor time_window_minutes ...
        sample_window Hz file
    
    fileName = files(file).name
    load(fileName)
    
    st = dm_spk_data.spike_times;
    % exclude overlapping templates
    
    numCells = length(st);
    spikeMatrix = zeros(sample_window,numCells);
    
    for cell = 1 : length(st)
        spikeTimes = round( st{cell} * 1000 , 0); % convert from s to ms with no d.p.
        spikeTimes = spikeTimes(spikeTimes<sample_window);
        %
        %     sum(spikeMatrix)
        spikeMatrix (spikeTimes,cell) = 1;
        progressbar(file / length(files) , cell / length(st))
    end
    FR_per_cell =  sum(spikeMatrix) ./ (length(spikeMatrix) / downHz) ;
    mean_FR_per_cell =  ( sum(sum(spikeMatrix)) / length(st) ) / (length(spikeMatrix) / downHz) ;
    lag_s = 0.05;
    
    tic
    adjM = getAdjM2(spikeMatrix, 'tileCoef', 0, lag_s, downHz);
    toc
    
    spikeMatrix = sparse(spikeMatrix);
    
    fs = downHz;
    xy = dm_spk_data.coordinates;
    
    save(strcat('MS_',fileName(1:end-4),'_revised_adjM.mat'),'st','xy','fs','adjM','spikeMatrix','-v7.3')
    
    disp('file done')
    progressbar(file / length(files) , cell / length(st))
end

clearvars -except files downHz downFactor time_window_minutes ...
    sample_window Hz file
%% plot adjMs
f1 = figure;

for file = 1 : length(files)
    fileName = files(file).name
    load(  strcat('MS_',fileName(1:end-4),'_revised_adjM.mat')  )

    subplot ( ceil( length(files) / 2 ) , 2 , file)
    imagesc(adjM)
    
    if file < length(files)
        clear adjM
    end
end

% add colorbar based on last culture
subplot ( ceil( length(files) / 2 ) , 2 , file +1)
addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'))
imagesc(adjM,'Visible','off'); aesthetics; axis off; 
cb = colorbar; cb.Location = 'Westoutside'; 
caxis( [min(min(adjM)),max(max(adjM))] ); cb.Label.String = 'Correlation';



%% Get distance matrix for each culture and put all in one array
% addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'))
% cd(dataDir);
files = dir('14*R.mat');

for file = 1 : length(files)
%     clearvars -except files file HD_MEA_Euclidean MEA_HD_binary
    fileName = files(file).name
    load(  strcat('MS_',fileName(1:end-4),'_revised_adjM.mat')  )
    %     remove NaNs etc
    adjM = weight_conversion(adjM,'autofix');
    
    %     threshold and binarise
    % adjM = threshold_proportional(abs(adjM), 0.2);
    
    % also remove <0.1 STTC 
    adjM = threshold_absolute(abs(adjM), 0.1);
    
    %     binarise
    adjM = weight_conversion(abs(adjM),'binarize');

    %     add network to cell array
    HD_MEA_binary{file}     = adjM;
    %     get distance mat and add to cell array
    distance_matrix         = squareform(pdist(xy,'euclidean')) ;
    HD_MEA_Euclidean{file}  = distance_matrix
    %     save cell array of adjMs and distance matrices
        
    clear adjM distance_matrix
end

save('HD_MEA_Euclidean_revised.mat','HD_MEA_Euclidean','-v7.3')
save('HD_MEA_binary_revised.mat','HD_MEA_binary','-v7.3')

% HD_MEA_binary_single(:,:,1)     = HD_MEA_binary{1};
% HD_MEA_Euclidean_single(:,:,1)  = HD_MEA_Euclidean{1};

% save('HD_MEA_binary_single.mat','HD_MEA_binary_single','-v7.3')
% save('HD_MEA_Euclidean_single.mat','HD_MEA_Euclidean_single','-v7.3')

%% plot distance matrices

load('HD_MEA_Euclidean_revised.mat')
f1 = figure;

for file = 1 : length(HD_MEA_Euclidean)
    d = HD_MEA_Euclidean{file};
    subplot ( ceil( length(HD_MEA_Euclidean) / 2 ) , 2 , file)
    imagesc(d)
    clear d
end

% add colorbar based on last culture
d = HD_MEA_Euclidean{file};
subplot ( ceil( length(HD_MEA_Euclidean) / 2 ) , 2 , file +1)
addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'))
imagesc(d,'Visible','off'); aesthetics; axis off; 
cb = colorbar; cb.Location = 'Westoutside'; 
caxis( [min(min(d)),max(max(d))] ); cb.Label.String = 'Distance \mum';

%% plot thresholded matrices
f1 = figure;
load('HD_MEA_binary.mat')
for file = 1 : length(HD_MEA_binary)
    subplot ( ceil( length(HD_MEA_binary) / 2 ) , 2 , file)
    imagesc(HD_MEA_binary{file})
    
    if file < length(HD_MEA_binary)
        clear adjM
    end
end

% add colorbar based on last culture
subplot ( ceil( length(HD_MEA_binary) / 2 ) , 2 , file +1)
addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'))
imagesc(HD_MEA_binary{end},'Visible','off'); aesthetics; axis off; 
c = parula;
colormap(c([1 256], :))
cb = colorbar; cb.Location = 'Westoutside'; 
caxis( [0 1] ); cb.Label.String = 'Correlation';
cb.Ticks = [0 1];
%% Run GM on HD MEA data
% addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output\GNM_Data'))
% load('MEAobserved28.mat')
% Run_MEA_generative_model
files = dir('14*.mat');

for file = 1 : length(files)
    fileName = files(file).name
    load(  strcat('MS_',fileName(1:end-4),'_spikes_adjM.mat')  )
%     remove NaNs etc
    adjM = weight_conversion(adjM,'autofix');
    original_adjM = adjM;
    
%     threshold and binarise
    adjM = threshold_proportional(abs(adjM), 0.2);
    thresholded_adjM = adjM;
    
%     binarise
    adjM = weight_conversion(adjM,'binarize');
%     figure;    imagesc(adjM)

% below measures all calculate in Run_MEAA_generative_model.m
%     distance_matrix =           triu( squareform(pdist(xy,'euclidean')) .* adjM );
% %     plot correlation between distance and correlation weight
% %     figure; scatter( distance_matrix(distance_matrix~=0) , thresholded_adjM(distance_matrix~=0) , '.k' );
%     
%     HD_MEA_Data{file,1}{1,1} = degrees_und(adjM)';
%     HD_MEA_Data{file,1}{2,1} = clustering_coef_bu(adjM);
%     HD_MEA_Data{file,1}{3,1} = betweenness_bin(adjM)';
%     HD_MEA_Data{file,1}{4,1} = distance_matrix(distance_matrix~=0);



    clear adjM original_adjM
end










%%
fileName = '1449.mat';
load(strcat('MS_',fileName(1:end-4),'_spikes_adjM.mat'))

thrOption = 'absolute';
thr = '0.1';
type = 'weighted';
plot_thr_network_fcn_V3(thrOption,thr,type,xy)

location.x = xy(:,1);
location.y = xy(:,2);

addpath(genpath('extended GLM\'));
addpath(genpath('learn_basis\'));
% load('simulated data_50 neurons.mat') % load simulated data
% spikes = data.spk;
load(strcat('C:\Users\alexd\Downloads\',fileName))
load(strcat('MS_',fileName(1:end-4),'_spikes_adjM.mat'))
clear all
load(strcat('MS_',fileName(1:end-4),'_spikes_adjM.mat'))
fileName = '1449.mat';
load(strcat('MS_',fileName(1:end-4),'_spikes_adjM.mat'))
spikes = data.spk; %spike train
load('simulated data_50 neurons.mat') % load simulated data
spikes = data.spk; %spike train
location.x = data.xx; % neuron location (x coordinate)
location.y = data.yy; % neuron location (y coordinate)
location.x = xy(:,1);
location.y = xy(:,2);
a = triu(adjM);b=triu(distance);scatter(a(:),b(:))
a = triu(adjM); b=triu(distance); scatter(a(:),b(:))
a = triu(adjM); b=triu(distance); scatter(b(:),a(:))
X = learning_basis(CCG,ignore);
model_fits(pre) = extendedGLM(CCG(pre,:), X, distance(pre,:),hyperparameter,ignore(pre,:));
results = detect_cnx(model_fits,ignore,threshold);
hyperparameter.binsize = 5; % binsize of the CCG (ms)
hyperparameter.interval = 100; % interval of the CCG (ms) eg. interval = 50 means the interval of the correlogram is [-25,25] ms
hyperparameter.tau0 = 1; %ms
hyperparameter.eta_w = 5;
hyperparameter.eta_tau = 20;
hyperparameter.eta_dt_coeff = 2;
[CCG, t, distance, ignore] = generate_correlogram(spikes,sr,location,hyperparameter,ignore_index);
X = learning_basis(CCG,ignore);
NN = size(CCG,1);
for pre = 1:NN % use parfor for parallel computing
    fprintf('neuron %i ',pre)
    model_fits(pre) = extendedGLM(CCG(pre,:), X, distance(pre,:),hyperparameter,ignore(pre,:));
end
[CCG, t, distance, ignore] = generate_correlogram(spikes,sr,location,hyperparameter,ignore_index);
X = learning_basis(CCG,ignore);
NN = size(CCG,1);
pre = 1
%-- 18/11/2020 15:43 --%
lastwarn
[msg,warnID] = lastwarn
warning('off','stats:glmfit:IterationLimit');
%-- 18/11/2020 15:53 --%
warning(state)