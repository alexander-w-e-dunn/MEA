function [CI,adjMm,adjMci,shuffled_spikeMat_sparse_all,adjMi] = AdjMProbabilisticThreshMEA(events,adjM,repNum,randMethod,adjMmethod,downSample,lag)

% Author RCFeord, June 2020
%
% this function generates a synthetic matrix of binary data 
%
% synthetic data is created using original data and shuffling events for
% each cell using one of two methods: 'shuffle' or 'circular'. All
% synthetic data is created from existing events
%
% it then averages the functional connectivity established across n number
% of iterated synthetic matrices

% Inputs:
%   events = binary matrix of spikes/neuronal events, columns are nodes/cells,
%            rows correspond to time points
%   adjM = adjacency matrix of real data
%   repNum = number of repetitions for the generation of synthetic datasets
%   randMethod = method used to create synthetic data: 'shuffle' 
%                (randomised) or 'circular' (conveyor belt)
%   adjMmethod = correlation method for adjacency matrix generation: 
%            'tileCoef' = STTC, 'correlation','partialcorr','xcorr'
%   downSample = downsampling (if none = 0)
%   lag = lag for the correlation or STTC

% Outputs:
%   CI = confidence interval for each edge on the basis of the
%        probabilistic edge weigth from 100 iterations of syntetic data
%   adjMm = mean edge weigth from 100 iterations of syntetic data              
%   adjMci = real adjacency matric thresholded at 95 confidence interval of
%            probabilistic edge weight

% updates AWE Dunn Oct 2020:
% output the randomised adjMs and the shuffled spike trains (all
% iterations). 
% clips from random point in time rather than random event
% 

for n = 1:repNum
    
    SynthDatBin = zeros(size(events));
    tic
    for i = 1:size(events,2)

        if strcmp(randMethod, 'shuffle')
            BA = find(events==1); % find burst activity
            BAi = BA(2:end)-BA(1:end-1);
            BAs(1) = BA(1); BAs = [BAs; BA(find(BAi>1)+1)]; % find end of burst
            BAe = BA(BAi>1); BAe = [BAe; BA(length(BA))]; % find end of burst
            BAl = BAe-BAs; % length of each burst
            BAlrand = BAl(randperm(length(BAl)));
            locs = randi([1 size(events,1)-20],size(BAl,2),1);
            for p = 1:length(locs)
                SynthDatBin(i,locs(p):locs(p)+BAlrand(p)) = 1;
            end
            clear BAi BAs BAe BAl BAlrand locs
        end
        
        if strcmp(randMethod, 'circular')

            if n == 1
                BAt = find(events(:,i)==1); % find burst activity
                VN = strcat('BA.i',num2str(i));
                eval([VN '= BAt;']);
                clear BAt
            end
            eval(['BAx = BA.i' num2str(i) ';']);
            if isempty(BAx)
                continue
            else
%                 locs = BAx(randi(numel(BAx)));
                locs = randi(length(events));
                SynthDatBin(1:size(events,1)-locs+1,i) = events(locs:size(events,1),i);
                SynthDatBin(size(events,1)-locs+2:size(events,1),i) = events(1:locs-1,i);
            end
            
            clear BAx locs
        end
        
    end
    toc
    
    % generate adjacency matrix from synthetic data
%     tic
    adjMs = getAdjM(SynthDatBin, adjMmethod, downSample, lag);
%     toc
    
    adjMs(adjMs<0) = 0;
    adjMs(isnan(adjMs)) = 0;
    for i = 1:length(adjMs)
        adjMs(i,i) = 0;
    end
    
    % save the synth adjMs 
    adjMi(:,:,n) = adjMs;
    % save the synth spike mat
    shuffled_spikeMat_sparse_all{n,1} = sparse(SynthDatBin);


end

% plot distribution of synthetic edge weights at 5 edges in the to 50% of
% non-zero connections from real matrix
adjMord = adjM(:);
adjMord = sort(adjMord);
adjMord(adjMord==0) = [];
sz = length(adjMord);
adjMord(1:sz/2) = [];
adjMord(isnan(adjMord)) = [];
randEdges = randi(length(adjMord),[5,1]);
randEdges = adjMord(randEdges);

p = [50 150 1200 300];
set(0, 'DefaultFigurePosition', p)
F1 = figure;
for p = 1:5
    subplot(1,5,p)
    [r,c] = find(adjM==randEdges(p));
    dat = adjMi(r(1),c(1),:);
    dat = dat(:);
    xbins = 0:0.05:1;
    hist(dat,xbins)
    hold on 
    plot([adjM(r(1),c(1)) adjM(r(1),c(1))],[0 1000],'r')
    if p == 1
    title('distribution of synthetic edge weights (red line real data)')
    end
    xlabel('edge weight')
    ylabel('frequency')
    xlim([0 1])
    ylim([0 100])
end

% average all synthetic adjacency matrices
adjMm = mean(adjMi,3);

% mean and standard dev per edge
meanCorrVal = mean(adjMm(:));

% confidence interval
CI = [];
adjMci = adjM;
% threshold the real adjM based on the 95% CI of the synthetic data
for i = 1:length(adjMm)
    for j = 1:length(adjMm)
        stdCorrVal = std(adjMi(i,j,:));
        SEMCorrVal = stdCorrVal/sqrt(repNum);
        ts = tinv([0.025 0.975], repNum-1);
        CI(i,j,[1 2]) = meanCorrVal + ts*SEMCorrVal;
        if CI(i,j,2)>adjM(i,j) 
            adjMci(i,j) = 0;
        end
    end
end

% plot adjacency matrix pre and post thresholding
p = [50 150 900 300];
set(0, 'DefaultFigurePosition', p)
F2 = figure;
subplot(1,2,1)
imagesc(adjM)
title('original adjacency matrix')
xlabel('nodes')
ylabel('nodes')
c = colorbar;
c.Label.String = 'correlation coefficient';
subplot(1,2,2)
imagesc(adjMci)
title('thresholded adjacency matrix')
xlabel('nodes')
ylabel('nodes')
c = colorbar;
c.Label.String = 'correlation coefficient';


end
