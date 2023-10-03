clear all

% load data
data=load('C:\Users\alexd\Downloads\170612c1fr4_01May2022.mat')

% set up filter and filter
highpass = 0.01 % in Hz
filterOrder = 3
wn = highpass / (data.Params.fs / 2);
[b, a] = butter(filterOrder, wn,'high');
filtdata = filtfilt(b,a,data.RawFluoro(:,1));

% plot raw
tiledlayout(3,1)
nexttile
plot(data.RawFluoro(:,1))
title('Raw')
% plot filtered
nexttile
plot(filtdata)
title('Filtered')
linkaxes
% plot filtered data normalised between 0 and 1
nexttile
norm_data = (filtdata - min(filtdata)) / ( max(filtdata) - min(filtdata) );
plot(norm_data)
title('Normalised')

