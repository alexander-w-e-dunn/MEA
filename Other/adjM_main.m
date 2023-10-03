%% get adjacency matrices
close all; clearvars -except refPeriod_ms

files = dir('*slice*mSpikes_3.mat');                   % where your .mat files are
files = files(~contains({files.name}, 'TTX'));
files = files(~contains({files.name}, 'ttx'));
files = files(~contains({files.name}, 'stim'));
files = files(~~contains({files.name}, 'FTD'));
files = files(~contains({files.name}, 'adjM'));
files = files(6:end);
%%%%%%% set these parameters first:
method = 'tileCoef';
sync_win_s = 0.200; %synchroncity window in seconds; e.g. 1 is +/- 1 s is considered synchronous (2DeltaT = 2 s)
rec_length_s = 360;
fs = 25000;
rec_length_samps = fs * rec_length_s;

%%%% downsampling:
num_samples = 1000 * rec_length_s; %defaults to 1 sample per ms
ds_rate = num_samples/rec_length_samps; %scalar for time_window %downsampling rate
sync_win = sync_win_s * ds_rate; %downsample the sync_win


fprintf(strcat('\n','\n','getting adjacency matrices...','\n','\n'))

% batch_getAdj_fcn(method,files,sync_win,num_samples,ds_rate);
% if recording lengths are all the same, uncomment above line and set
% rel_length_s above; otherwise, use below code to automatically find
% recording lengths.

for i = 1:length(files)
    files1 = files(i);
    spikes = struct2cell(load(files(i).name,'*Spikes'));
    rec_length_s = length(spikes{1}) / fs;
    rec_length_samps = fs * rec_length_s;
    num_samples = 1000 * rec_length_s; %defaults to 1 sample per ms
    ds_rate = num_samples/rec_length_samps; %scalar for time_window %downsampling rate
    sync_win = sync_win_s * ds_rate; %downsample the sync_win
    
    batch_getAdj_fcn(method,files1,sync_win,num_samples,ds_rate);
end