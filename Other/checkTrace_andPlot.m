% plotting activity to exclude

%% load
load('MPT190403_2A_DIV28.mat')
%load spikes
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\cSpikes.2507'
load('MPT190403_2A_DIV28_cSpikes_L0.2507.mat')
%load filtmat
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\filt_mats'
load('MPT070119_2A_DIV28_Filtd.mat')

data=zeros(size(dat));
data=dat-mean(dat); %normalise


ch71=dat(:,find(channels==71));
plot(ch71-mean(ch71))

ch72=dat(:,find(channels==72));
plot(ch72-mean(ch72))

ch82=dat(:,find(channels==82));
plot(ch82-mean(ch82))


%% grid trace plot
%currently doesnt work with raw data; adjust scales
downFactor = 1000; % down sample factor for making grid trace
cd 'D:\MECP2_2019_AD\Scripts_and_Output'
%filteredMatrix=filteredMatrix(timeRange,:);

gridTrace_AD(data(:,:), downFactor,'id',channels,fs) %need to add folders to path containing aesthetics function etc.

%% mutliple single trace plots
yGap = 100; % vertical gap bewteen traces 
electrodesToPlot = [51:60]; % list of electrodes to plot
timeRange = 1: fs * 720;
%faster if downsample:
traces = downsample(dat, downFactor);
timeRange = timeRange(1:end/1000);

figure 

for electrode = 1:length(electrodesToPlot)
            
            plot(traces(timeRange, electrodesToPlot(electrode)) - yGap * (electrode -1)...
                ,'color','b')
            hold on 
           
end

aesthetics 
removeAxis 
scalebar

%% plot 4 plots

%ensure from top, raw data loaded and normalised, filt data loaded, spike
%data loaded

%data=dat-mean(dat); %ensure data normalised

chanID = 66;
chan = find(channels==chanID);

figure
fs=25000;

% 10s raw

timeRange = 1:10*fs;
trace = data(timeRange,chan);

subplot(4,1,1)
plot([1:length(trace)]/fs,trace)
aesthetics
removeAxis 

       sb = scalebar;
       %sb.YLen = round(max(trace)/2); %round to nearest 10
       sb.YLen = 20; %round to nearest 10
       sb.XLen = 1; 
       sb.YUnit = '\muV';
       sb.XUnit = 's'; 
       sb.Position = [-0.1, -50];
       %may need to adjust accoridng to legnth of plot window
       sb.hTextX_Pos = [0.4, -10];
       sb.hTextY_Pos = [-0.2, 5];
       
       
% 10s filtered

filttrace = filteredData(timeRange,chan);

subplot(4,1,2)
plot([1:length(filttrace)]/fs,filttrace)
aesthetics
removeAxis 

       sb = scalebar;
       %sb.YLen = round(max(trace)/2); %round to nearest 10
       sb.YLen = 20; %round to nearest 10
       sb.XLen = 1; 
       sb.YUnit = 'a.u.';
       sb.XUnit = 's'; 
       sb.Position = [-0.1, -15];
       %may need to adjust accoridng to legnth of plot window
       sb.hTextX_Pos = [0.4, -10];
       sb.hTextY_Pos = [-0.2, 5];
       linkaxes

       
       %% sep figure 10ms
       figure
% 40ms during spikes raw dat
spikes = full(cSpikes);
spikeT = spikes(:,chan);
%examine spikeT to get good set of spikes if spike count low
spikeTrace = data(614001:615000,chan);

subplot(4,1,3)
plot([1:length(spikeTrace)]/(fs/1000),spikeTrace)
aesthetics
removeAxis 

       sb = scalebar;
       %sb.YLen = round(max(trace)/2); %round to nearest 10
       sb.YLen = 20; %round to nearest 10
       sb.XLen = 1; 
       sb.YUnit = '\muV';
       sb.XUnit = 'ms'; 
       sb.Position = [-0.1, -30];
       %may need to adjust accoridng to legnth of plot window
       sb.hTextX_Pos = [0.4, -5];
       sb.hTextY_Pos = [-1, 5];
       

% 40ms during spikes filtered

spikeTraceF = filteredData(614001:615000,chan);

subplot(4,1,4)
plot([1:length(spikeTraceF)]/(fs/1000),spikeTraceF)
aesthetics
removeAxis 

       sb = scalebar;
       %sb.YLen = round(max(trace)/2); %round to nearest 10
       sb.YLen = 10; %round to nearest 10
       sb.XLen = 1; 
       sb.YUnit = 'a.u.';
       sb.XUnit = 'ms'; 
       sb.Position = [-0.1, -11];
       %may need to adjust accoridng to legnth of plot window
       sb.hTextX_Pos = [0.4, -5];
       sb.hTextY_Pos = [-1, -1];
       
       linkaxes
       %ylim([-10,10]);
       
%add a point to indicate the spike time
sTime = (614976-614001)/(fs/1000);
x = sTime;
y = -15;
plot(x,y,'^','LineWidth',5,'MarkerFaceColor','r','MarkerEdgeColor','r');

%% add in plot channel with spike points marked underneath
option = 'line' %line or dot
scb = 1; %0 for no scalebar 1 for include scalebar

%channelIDofInterest= 27;
%chan_num = find(channels==27);
chan_num = 27;
timetoplots = 1; %in seconds
fs=25000;
timewindow = 1:timetoplots*fs;

spikeMat = full(cSpikes);
spikeTrain1 = spikeMat(:,chan_num);

figure
%plot raw data
subplot(4,1,1)
raw_trace=dat(timewindow,chan_num)-mean(dat(:,chan_num));
x_time_s=[1:length(timewindow)]/(fs);
plot(x_time_s',raw_trace);
aesthetics
axis off
if scb == 1
sb1 = scalebar
       sb1.YLen = 20; %round to nearest 10
       %sb.XLen = length(electrodeMatrix)/fs; 
       sb1.YUnit = '\muV';
       %sb.XUnit = 'ms';
       sb1.XUnit = 's';
       sb1.Position = [-timetoplots/20, min(raw_trace)-7];
       %may need to adjust accoridng to legnth of plot window
       sb1.hTextX_Pos = [timetoplots/20,-sb1.YLen/2];
else
end
 title(strcat(int2str(timetoplots),'(s) of raw voltage trace'))
       
       
%plot filtered data


subplot(4,1,2)

%get MS spikes -5SD & filt dat
method = 'Manuel'; multiplier = 5; L=0;
[spikeTrain, finalData, threshold] = detectSpikes(dat(timewindow,chan_num), method, multiplier, 0);
x_time_s=[1:length(timewindow)]/(fs);
plot(x_time_s',finalData);
aesthetics
axis off
title(strcat(int2str(timetoplots),'(s) of filtered voltage trace'))
linkaxes

if scb == 1
sb2 = scalebar
       sb2.YLen = 20; %round to nearest 10
       %sb.XLen = length(electrodeMatrix)/fs; 
       sb2.YUnit = 'a.u.';
       %sb.XUnit = 'ms';
       sb2.XUnit = 's';
       sb2.Position = [-timetoplots/20, min(raw_trace)-7];
       %may need to adjust accoridng to legnth of plot window
       sb2.hTextX_Pos = [timetoplots/20,-sb1.YLen/2];
       
else
end
%add threshold line
threshold_trace = ones(length(timewindow),1)*threshold;
hold on
plot(x_time_s',threshold_trace,'color','r');
hold off
        
%plot train method one
subplot(4,1,3)
singleRastPlot(spikeTrain1(timewindow,1), option) 
aesthetics
title('Wavelet-based method')
%sb3 = scalebar
       %sb3.Position = [-timetoplots/20, min(raw_trace)-7];

%plot spike train method 2
subplot(4,1,4)
singleRastPlot(spikeTrain(timewindow,1), option) 
aesthetics
title('Schroeter et al (2015) method')
%sb4 = scalebar
       %sb4.Position = [-timetoplots/20, min(raw_trace)-7];






