clear all
% %  example
% clearvars; clc;
% dataPath = '/Users/jeremi/mea/data/organoid/';
% savePath = '/Users/jeremi/mea/data/organoid/';
% option = 'list';
% files = {'200617_slice1.mat'};
% load params
% params.costList = -0.4;
% % params.subsample_time = [1, 60];
% 
% batchDetectSpikes(dataPath, savePath, option, files, params);
%%
% custom_spike_scriptsDir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts\';
custom_spike_scriptsDir = 'C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts\Other_scripts\WATERS-master\';
addpath(genpath(custom_spike_scriptsDir));
% dataPath = ['C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\']; % – set path to folder with data
% dataPath = ['D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\By_culture\MPT\Recordings\']; % – set path to folder with data
dataPath = ['D:\MECP2_2019_AD\Scripts_and_Output\2.File_Conversion_Output'];
savePath = ['C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Output\Mecp2\']; % – set path to output folder
cd(custom_spike_scriptsDir)

% setParams(); %– set spike detection parameters

% input cell array containing file name strings
% filestruct = dir([dataPath  '\' 'MPT200205_4A*.mat']); files = struct2cell(filestruct)'; files = files(:,1);
filestruct = dir([dataPath  '\' '1907*.mat']); files = struct2cell(filestruct)'; files = files(:,1);

% 'list' option for cell array, 'path' option for all files in a directory
option = 'list';
load params

cd(dataPath)

tic
% getSpikesTS(dataPath, savePath);
% batchDetectSpikes(dataPath, savePath, files);
batchDetectSpikes(dataPath, savePath, option, files, params);
toc

%% check if spikes found
% cd(savePath)
% find files in save path
files = dir([savePath  '\' '*spikes.mat*']);
files = files(~contains({files.name}, 'Spikes','IgnoreCase',false));
% files = files(1);
% sort files according to date and get latest file
[~,index] = sortrows({files.date}.'); files = files(index(end:-1:1)); clear index
load(files(1).name);

for elec = 1 : length(spikeTimes)
    % counts(elec) = sum(spikeTimes{elec}.mea) ;
    % counts(elec) = length(spikeTimes{elec}.mea) ;
    counts(elec) = length(spikeTimes{elec}.bior1p5) ;
end

sum(counts) 

MEAgraphics(1)
histfit(log10(counts(counts>0)),20,'kernel');xlim([0 4]);xticks(0:4);xticklabels(10.^xticks);aesthetics
xlabel('Spike count'); ylabel('Frequency (electrodes)');

%%
% spikeMatrix