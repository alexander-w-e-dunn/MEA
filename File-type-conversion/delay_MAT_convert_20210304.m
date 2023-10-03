% set number of hours to delay onset of raw to mat conversion
hours         = 3; 
% directories
ext           = '.raw';
% scriptsFolder = 'C:\Users\Windows\Documents\MATLAB';
scriptsFolder = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';

% inputFolder   = 'G:\Alex\RAW files\Collaboration_with_Paul';
% outputFolder  = 'F:\Alex\MAT files\Collaboration_with_Paul';

% inputFolder   = 'C:\Users\Windows\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Mecp2 Project Jan 2019-\RAW files';
% outputFolder  = 'C:\Users\Windows\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Mecp2 Project Jan 2019-\MAT files';

% inputFolder   = 'C:\Users\Windows\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Nanobodies';
% outputFolder  = 'C:\Users\Windows\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Nanobodies';

inputFolder   =  'C:\Users\alexd\Downloads';
outputFolder  =  'C:\Users\alexd\Downloads';

cd(scriptsFolder);
addpath(genpath(scriptsFolder));
%% create timer
% input #2 to below timer function is the delay in s
T = timer('StartDelay',hours*60*60,'TimerFcn',...
    @(src,evt)MEAbatchConvert_20210304(ext,inputFolder,outputFolder,scriptsFolder)); 
% run
start(T)

%% create timer
% % input #2 to below timer function is the delay in s
% T2 = timer('StartDelay',hours*2*60*60,'TimerFcn',...
%     @(src,evt)batchGetSpike); 
% % run
% start(T2)