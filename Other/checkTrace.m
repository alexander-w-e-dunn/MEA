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
electrodesToPlot = [15,21,32]; % list of electrodes to plot
timeRange = 1: fs * 720;
%faster if downsample:
traces = downsample(data, downFactor);
timeRange = timeRange(1:end/1000);

figure 

for electrode = 1:length(electrodesToPlot)
            
            plot(traces(timeRange, electrodesToPlot(electrode)) - yGap * (electrode -1))
            hold on 
           
end

aesthetics 
removeAxis 

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

