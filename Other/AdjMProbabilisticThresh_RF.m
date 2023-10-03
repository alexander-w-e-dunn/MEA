function [adjMci] = AdjMProbabilisticThresh_v3(events,adjM,repNum,adjMmethod,downSample,lag,fs,ploton,tail)

% Updated by HSmith, Cambridge, December 2020
% Author RCFeord
%
% This function generates a synthetic matrix of binary data
%
% Synthetic data is created using original data and shuffling events for
% each cell using 'circular' method. All synthetic data is created from 
% existing events
%
% It then averages the functional connectivity established across n number
% of iterated synthetic matrices
%
% INPUTS:
%   events = binary matrix of spikes/neuronal events, columns are nodes/cells,
%            rows correspond to time points
%   adjM = adjacency matrix of real data
%   repNum = number of repetitions for the generation of synthetic datasets
%   adjMmethod = correlation method for adjacency matrix generation: 
%            'tileCoef' = STTC, 'correlation','partialcorr','xcorr'
%   downSample = downsampling (if none = 0)
%   lag = lag for the correlation or STTC
%   fs = sampling frequency
%   ploton:
%       = 0 do not show figures
%       = 1 show figures
%   tail = confidence interval. Eg input 0.05 for p = 0.05 thresholding and
%          0.025 for p = 0.025 thresholding. No default set


% OUTPUTS:         
%   adjMci = real adjacency matrix thresholded at "tail" confidence interval
%            of probabilistic edge weight

%% Generate synthetic data
num_frames = size(events,1);
num_nodes = size(events,2);

for n = 1:repNum
    
    % Create a matrix the same size as the real data matrix ('events')
    SynthDatBin = zeros(size(events));
    
    locs = randi(num_frames,1,num_nodes);
       
    for i = 1:num_nodes
        
        SynthDatBin(1 : end - locs(i) +1 , i) = events(locs(i) : end , i);
        SynthDatBin(end - locs(i) +2 : end , i) = events(1 : locs(i) -1 , i);
        
    end
        
    % Generate adjacency matrix from synthetic data
    adjMs = getAdjM_CaPop(SynthDatBin, adjMmethod, downSample, lag, fs);

    % Add to stack of synthetic adjacency matrices
    adjMi(:,:,n) = adjMs; 
end

% Remove negatives, NaNs, and nonzero diagonals
adjMi(adjMi<0) = 0;
adjMi(isnan(adjMi)) = 0;
adjMi = bsxfun(@times, ~eye(num_nodes), adjMi);

%% Plot some synthetic data to validate
% Plot distribution of synthetic edge weights at 5 edges in the top 50% of
% non-zero connections from real matrix
if ploton == 1
    adjMord = adjM(:); % Turn matrix into vector 
    adjMord = sort(adjMord); % Lowest to highest 
    adjMord(adjMord==0) = []; % Remove 0 values
    sz = length(adjMord);
    adjMord(1:sz/2) = []; % Take top 50% data
    adjMord(isnan(adjMord)) = []; % Remove NaNs
    randEdges = randi(length(adjMord),[5,1]); % Choose 5 random points
    randEdges = adjMord(randEdges); % and select their values in adjMord

    p = [50 150 1200 300];
    set(0, 'DefaultFigurePosition', p)
    F1 = figure;
    for p = 1:5
        subplot(1,5,p)
        [r,c] = find(adjM==randEdges(p)); % Locate one of the 5 random points in adjM 
        dat = adjMi(r(1),c(1),:); % Run along all repNum points  
        dat = dat(:); % Vectorise 
        xbins = 0:0.05:1;
        hist(dat,xbins) % Histogram frequency of occurrence of those random points
        hold on 
        plot([adjM(r(1),c(1)) adjM(r(1),c(1))],[0 1000],'r') % Plot value in adjM
        if p == 1
        title('Distribution of synthetic edge weights (red line real data)')
        end
        xlabel('edge weight')
        ylabel('frequency')
        xlim([0 1])
        ylim([0 100])
    end
end

%% Stats test         
% Threshold each element if >= top "tail" % of data 
adjMci = adjM;
cutoff_point = ceil((1 - tail) * repNum);
for i = 1:size(adjMi,1)
    for j = 1:size(adjMi,2)
        Eu = sort(adjMi(i,j,:),'ascend');
        if Eu(cutoff_point) > adjM(i,j)
            adjMci(i,j) = 0;
        end
    end
end

%% Plot adjacency matrix pre and post thresholding
if ploton == 1
    p = [50 150 900 300];
    set(0, 'DefaultFigurePosition', p)

    F2 = figure;

    subplot(1,3,1)
    imagesc(adjM)
    title('Original adjacency matrix')
    xlabel('nodes')
    ylabel('nodes')
    c = colorbar;
    c.Label.String = 'correlation coefficient';

    subplot(1,3,2)
    imagesc(adjMci)
    title('Thresholded adjacency matrix')
    xlabel('nodes')
    ylabel('nodes')
    c = colorbar;
    c.Label.String = 'correlation coefficient';

    subplot(1,3,3)
    imagesc(adjM - adjMci)
    title('Difference between original and thresholded')
    xlabel('nodes')
    ylabel('nodes')
    c = colorbar;
    c.Label.String = 'correlation coefficient';

    limits = F2.Children(5).Limits;
    F2.Children(3).Limits = limits;
    F2.Children(1).Limits = limits;
    F2.Children(2).CLim = limits;
    F2.Children(4).CLim = limits;
end

end