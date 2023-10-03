% Load data;
addpath(genpath('C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\Organoid\Spikes\'));
load('191210_slice1_DIV_g04_2018_L_-0.0627_spikes.mat');

% Some params and pre-allocation
fs = 25000;
lag = 10; % [ms]
start_time = 0; end_time = spikeDetectionResult.params.duration;

%%
spike_matrix = spikeTimeToMatrix(spikeTimes, start_time, end_time, fs);
%%
%downsample
fs = 1000;
adjM2 = getAdjM2(spikes,method,downSample,lag,fs)