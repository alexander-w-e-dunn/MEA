clear all

scriptsDir  = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts\';
% set data dir
dataDir     = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\Mecp2\Spikes';
% fileName = 'MPT200108_2B_DIV28_cSpikes_L-0.0627_RF1.mat';
fileName    = 'MPT190403_2B_DIV21_cSpikes_L-0.0627_RF1.mat';
% fileName = 'C:\Users\alexd\Downloads\1449.mat';

%% compute adjM
addpath(genpath(scriptsDir));
cd(dataDir);
load(fileName); 

spikeMatrix = full(cSpikes);
lag = 0.05;
fs = 25000;
parallel = 0;
% runs faster if you run on sparse matrix
adjM = sttc_fcn(cSpikes,lag,fs,parallel);
% figure; imagesc(abs(adjM)); colorbar; caxis([0 1]);

%% threshold matrix A) n permutations
numPermutations = 100;

new_sample_rate = 1000;
downFactor = fs/new_sample_rate;
downMatrix = downSampleSum(spikeMatrix, length(spikeMatrix)/downFactor);
downMatrix(find(downMatrix>1)) = 1;

%% load adjacency matrix and downsample
addpath(genpath(scriptsDir));
cd(dataDir);
load(fileName); load(strcat(fileName(1:end-4),'_adjM_0.05.mat'));
numPermutations = 100;
sync_win_s = 0.05;
spikeMatrix = zeros(size(cSpikes));
spikeMatrix = full(cSpikes);
% spikeTimes = find(spikeMatrix==1);
if ~exist('fs')
    fs = 25000;
end
new_sample_rate = 1000;
downFactor = fs/new_sample_rate;
downMatrix = downSampleSum(spikeMatrix, length(spikeMatrix)/downFactor);
clear spikeMatrix
downMatrix(find(downMatrix>1)) = 1;
%% run corelation control
tic
for permutation = 1 : numPermutations
    % tic
    synth_mat_i = zeros(size(downMatrix));
    for elec = 1 : size (downMatrix,2)
        randomSample = randi(length(downMatrix));
        
        synth_mat_i(1 : length(downMatrix) - randomSample , elec)...
            = downMatrix(randomSample+1 : size(downMatrix,1),elec);
        
        synth_mat_i(size(downMatrix,1) - randomSample + 1 : size(downMatrix,1),elec)...
            = downMatrix(1 : randomSample, elec);
        
        clear randomSample
    end
    % toc
    synth_mats{permutation} = synth_mat_i;
    clear synth_mat_i 
end
toc; sound(sin(1000:3000))

% v1 = sum(synth_mats{1});
% v2 = sum(downMatrix); v2-v1

% downSampleRate = new_sample_rate/fs;

%% correlate adjMs
tic
for permutation = 1 : numPermutations
%     tic
%     synth_adjMs{permutation} = getAdjM(synth_mats{permutation}, 'tileCoef', 0, sync_win_s * downSampleRate);
%     toc
    Sspikes = sparse(synth_mats{permutation});
    SadjM = sttc_fcn(Sspikes,lag,fs,parallel);
    synth_adjMs{permutation} = SadjM;
    clear Sspikes SadjM
end
toc;sound(sin(1000:3000))

%% save and analyse
save(strcat(fileName(1:end-4),'_ctrl.mat'),'synth_mats','synth_adjMs','channels','-v7.3')

matrix = synth_adjMs{1};
v = find(matrix > 0.1 & matrix < 0.9);
edge = v(1);

for i=1:length(synth_adjMs)
    matrix = synth_adjMs{i};
    edges(i) = matrix(edge);
end
figure; histfit(edges);%sound(sin(1000:3000))
adjM(edge)
matrix = synth_adjMs{1};matrix(edge)