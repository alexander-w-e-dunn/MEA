clear all
load('200127_FTDOrg_GrpA_3A_Slice8_mSpikes_3.mat')
load('200127_FTDOrg_GrpA_3A_Slice8_mSpikes_3_adjM.mat');

%%
events = full(mSpikes);
% events = full(mSpikes)';
fs=25000;
rec_duration = 360; % in seconds
events = events(1:rec_duration*fs,:);
% repNum = 100;
repNum = 3;
randMethod = 'circular';
adjMmethod = 'tileCoef';
rec_s = length(events)/fs;
downSample = rec_s*1000; % new number of sample
ds_rate=downSample/length(events);
lag_s = 0.175;
lag = ds_rate * lag_s;

if downSample ~= length(events)
    events_ds = downSampleSum(events,downSample);
end
adjM1 = getAdjM(events_ds, adjMmethod, downSample, lag);
fpos = [50 150 1200 300];
set(0, 'DefaultFigurePosition', fpos);
% confirmed that this produces the expected adjM!
% f1=figure; subplot(1,2,1); imagesc(adjM); subplot(1,2,2); imagesc(adjM1)

%% create synth mat
% create vector of random times
% r_vec = rand(1,60) * length(events_ds);
rand_events = zeros(size(events_ds));
for elec = 1:size(events,2)
   rand_events(1 : r_vec(elec) , elec) = events(r_vec(elec) : end , elec);
   rand_events(r_vec(elec)+1 : end , elec) = events(1 : r_vec(elec)-1 , elec);
end