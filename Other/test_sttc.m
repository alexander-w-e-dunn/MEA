clc
%% IMPORTANT: if you don't have the appropriate .mex version of your .c 
%             script (e.g. .mexmaci64 for OSX), you need to first run (in 
%             command window):
%             mex sttc.c -R2018a

% @author: JJChabros, 2021

% mex sttc_jeremi.c -R2018a

% Load data;
addpath(genpath('C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\Organoid\Spikes\'));
load('191210_slice1_DIV_g04_2018_L_-0.0627_spikes.mat');

% Initialize parallel computing cluster
% JJC: Dunno if it'll make any difference to the execution time but let's
% give it a go!
poolobj = gcp('nocreate');
delete(poolobj);
parpool('local', 4); % number of cores

% Some params and pre-allocation
fs = 25000;
lag = 0.01; % [s]
numChannel = length(channels);
combChannel = nchoosek(1:numChannel, 2);
A = zeros(1, length(combChannel));
adjM = NaN(numChannel, numChannel);
%% Run sttc

% NOTE: if you don't want to restart MATLAB and lose all of your progress,
% it is paramount that the arguments passed to the sttc function are of
% appropriate type (as specified below)

% TODO: make it a flexible-input function
tic
parfor i = 1:length(combChannel)
% for i = 1:length(combChannel)
    dtv = lag; % [s]
    spike_times_1 = double(spikeTimes{combChannel(i,1)}.bior1p5/25000);
    spike_times_2 = double(spikeTimes{combChannel(i,2)}.bior1p5/25000);
    N1v = int16(length(spike_times_1));
    N2v = int16(length(spike_times_2));
    dtv = double(dtv);
    Time = double([0 spikeDetectionResult.params.duration]);
    tileCoef = sttc_jeremi(N1v, N2v, dtv, Time, spike_times_1, spike_times_2);
%     coder.ceval('sttc_jeremi',N1v, N2v, dtv, Time, spike_times_1, spike_times_2);
    row = combChannel(i,1);
    col = combChannel(i,2);
    A(i) = tileCoef; % Faster to only get upper triangle so might as well store as vector
end
toc

% Vector -> matrix
for i = 1:length(combChannel)
    row = combChannel(i,1);
    col = combChannel(i,2);
    adjM(row, col) = A(i);
    adjM(col, row) = A(i);
end


