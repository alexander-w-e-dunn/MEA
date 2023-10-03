% INPUTS:
%     fileName        = spikematrix file (sparse(num samples x num channels))
%     lag             = sync. window in s
%     numPermutations = number of permutations for control deafult is 1000
%     method          = 'CI','PCT', or 'AUC' string input for confidence
%                         interval, percentile or area under curve threshold respectively,
%                         default is 'PCT'
%     parameter       = specifies value for thresholding, e.g.
%                         input 95 for 95th percentile or 95 % CI.
%
% OUTPUTS:
%     adjM            = the observed sttc matrix
%     thr_adjM        = thresholded matrix
%     thr             = vector of thresholds for each channel
%     synth_adjMs     = the synthetic adjMs 
%     channels        = vector of MEA channel IDs
%     EdgeIDs         = matrix giving each possible pair of channel IDs
%     EdgeNum         = matrix giving channel indeces for each edge
% 
% Authors AWE DUNN and R FEORD, U. Cambridge, 2021; 
% Uses functions by J CHABROS AND T SIT
% 
% EXAMPLE SETUP:
% 
% clear all
% % set directories
% scriptsDir  = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts\';
% dataDir     = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\Mecp2\Spikes';
% addpath(genpath(scriptsDir));
% cd(dataDir);
% 
% lag = 0.05;
% fs = 25000;
% numPermutations = 100;
% method = 'PCT';
% alpha = 90; % 95 means 5 % false positive rate / sig. at p0.05 level
%             % could use 90 as this is comparable to 90 % proportion
%             % threshold
% 
% TO DO:
% - add CI and area under curve options
% - make this work for other correlation metrics
% - adapt for directed matrices

function [adjM,thr_adjM,thr,synth_adjMs,EdgeIDs,EdgeNum] = threshold_MEA_STTC(sparse_spikeMat,lag,fs,numPermutations,method,alpha,channels)


%% compute adjM
% load(fileName); 
spikeMatrix = full(sparse_spikeMat);
% parallel computing? 0 = no, 1 = yes.
parallel = 0;
% runs faster if you run on sparse matrix
adjM = sttc_fcn(sparse_spikeMat,lag,fs,parallel);
% figure; imagesc(abs(adjM)); colorbar; caxis([0 1]);

%% downsample spike matrix 
% downsampled matrix used to generate synthetic data at 1 sample/ms
new_sample_rate = 1000;
downFactor = fs/new_sample_rate;
downMatrix = downSampleSum(spikeMatrix, length(spikeMatrix)/downFactor);
downMatrix(find(downMatrix>1)) = 1;
clear spikeMatrix
%% Create shuffled spike matrices
tic
for permutation = 1 : numPermutations
    disp(newline + "permutation: " + num2str(permutation) + newline)
    % tic
    synth_mat_i = zeros(size(downMatrix));
    for elec = 1 : size (downMatrix,2)
        % disp(newline + "electrode: " + num2str(elec) + newline)
        randomSample = randi(length(downMatrix));
        
        synth_mat_i(1 : length(downMatrix) - randomSample , elec)...
            = downMatrix(randomSample+1 : size(downMatrix,1),elec);
        
        synth_mat_i(size(downMatrix,1) - randomSample + 1 : size(downMatrix,1),elec)...
            = downMatrix(1 : randomSample, elec);
        
        clear randomSample
    end
    % toc
    synth_mats{permutation} = sparse(synth_mat_i);
    clear synth_mat_i 
end
disp(newline + "       --> spike matrix shuffles generated" + newline)
toc; % sound(sin(1000:3000))

%% correlate shuffled spike matrices
tic
for permutation = 1 : numPermutations
    disp(newline + "permutation: " + num2str(permutation) + newline)
%     tic
%     synth_adjMs{permutation} = getAdjM(synth_mats{permutation}, 'tileCoef', 0, sync_win_s * downSampleRate);
%     toc
    Sspikes = synth_mats{permutation};
    SadjM = sttc_fcn(Sspikes,lag,fs,parallel);
    synth_adjMs{permutation} = SadjM;
    clear Sspikes SadjM
end
disp(newline + "       --> shuffled spike matrices correlated" + newline)
toc;% sound(sin(1000:3000))

%% threshold the matrix and save shuffled spike matrics, corelation controls and thresholded matrix
% not saving synth. spike matrices due to file size
% to do: other thresholding options
% 1. CI 
% 2. top 5 % area under curve
for permutation = 1 : length(synth_adjMs)
    all_mats(:,:,permutation) = synth_adjMs{permutation};
end

thr_adjM = adjM;
EdgeIDs = nchoosek(channels,2);
EdgeNum = nchoosek(1:length(channels),2);

for edge          = 1 : length(EdgeNum)
    weights       = squeeze(all_mats(EdgeNum(edge,1),EdgeNum(edge,2),:));
    % get threshold by specified method
    if strcmp(method,'PCT')
        thr(edge) = prctile(weights,alpha);
        % add other thr. options here
    end
    % remove edges less than threshold
    if adjM(EdgeNum(edge,1),EdgeNum(edge,2)) < thr(edge)
        % threshold edge in upper triangle
        thr_adjM(EdgeNum(edge,1),EdgeNum(edge,2)) = 0;
        % threshold edge in lower triangle
        % to make this script work for a directed matrix, need to calculate
        % the threshold for thr_adjM(EdgeNum(edge,2),EdgeNum(edge,1)) i.e.
        % index 2, then index 1 which means re-running loop for all pairs
        % of EdgeIDs
        thr_adjM(EdgeNum(edge,2),EdgeNum(edge,1)) = 0;
    end
    clear weights
end

end