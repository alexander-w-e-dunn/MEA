% load spike data and adjacency matrix
clear all ; close all
dataDir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\AnalysisOutput\By_culture\MEC\Spikes\cSpikes_L-0.0627_RP1'
cd(dataDir)
scriptsDir      = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
addpath(genpath(scriptsDir));
load('MEC_220425_2A_DIV28_cSpikes_L-0.0627_RP1.mat')
load('MEC_220425_2A_DIV28_cSpikes_L-0.0627_RP1_adjM_0.02.mat')
adjM = weight_conversion(adjM,'autofix')
b = threshold_proportional(adjM,.1)
b = weight_conversion(b,'binarize')
nexttile
imagesc(adjM); caxis([0 1])
nexttile
imagesc(b)

%% plot the graph of a binary matrix
g=graph(b);
c = num2str(channels)
x = str2num(c(:,1))
y = str2num(c(:,2))
cs = string(channels)
% cm = brewermap(sum(b,'all')/2,'Greys')
cm = flip(gray(sum(b,'all')/2));
w = abs(adjM(triu(logical(b))));
d = degrees_und(b);
fr = full(sum(cSpikes))';
nexttile
p=plot(g,'XData',x,'YData',-y,'NodeLabel',cs,'EdgeColor',1-(w/max(w)*[1 1 1]),'LineWidth',2,'MarkerSize',d+1,'NodeColor',1-(fr/max(fr)*[0 1 1]))

aesthetics
axis off
