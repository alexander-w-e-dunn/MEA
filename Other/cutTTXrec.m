%% edit raw dat
clear all
fileName = '191210_slice1_DIV_g04_2018.mat'
load(fileName)
dat = dat(1:120*fs,:);
save(strcat(fileName(1:end-4), '_2min.mat'), 'ADCz','channels','dat','fs')
%% ttx edit raw dat
clear all
fileName = '191210_slice1_DIV_g04_2018_TTX.mat'
load(fileName)
dat = dat(1:120*fs,:);
save(strcat(fileName(1:end-4), '_2min.mat'), 'ADCz','channels','dat','fs')
%%
save(fileName, 'aSpikes','channels','thresholds','refPeriod_ms')