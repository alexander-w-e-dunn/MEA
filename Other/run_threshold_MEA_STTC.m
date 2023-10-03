clear all
% set directories
scriptsDir  = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts\';
dataDir     = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\Mecp2\Spikes';
addpath(genpath(scriptsDir));
cd(dataDir);

% fileName = 'MPT200108_2B_DIV28_cSpikes_L-0.0627_RF1.mat';
fileName    = 'MPT190403_2B_DIV21_cSpikes_L-0.0627_RF1.mat';
% fileName = 'C:\Users\alexd\Downloads\1449.mat';
lag = 0.05;
fs = 25000;
numPermutations = 100;
method = 'PCT';
alpha = 95; % 95 means 5 % false positive rate / sig. at p0.05 level
            % could use 90 as this is comparable to 90 % proportion
            % threshold
% run thresholding for each file 
tic
[adjM,thr_adjM,thr,synth_adjMs,channels,EdgeIDs,EdgeNum] = threshold_MEA_STTC(fileName,lag,fs,numPermutations,method,alpha);
toc
% figure; imagesc(adjM); colorbar; caxis([0 1]);
% figure; imagesc(thr_adjM); colorbar;caxis([0 1]);
% figure; imagesc(adjM-thr_adjM); colorbar;caxis([0 1]);
% diff = adjM-thr_adjM;
% nanmean(nanmean(diff(diff>0)))
% nanmean(nanmean(thr_adjM(thr_adjM>0)))
%%
save(strcat(fileName(1:end-4),'_adjM_thresholded.mat'),'adjM','thr_adjM','thr','synth_adjMs','channels','EdgeIDs','EdgeNum','-v7.3')
