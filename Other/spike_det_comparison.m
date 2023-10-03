% run spike_overlay_fcn to plot 30 ms plot for 'MPT200107_1A_DIV28_cSpikes_L-0.1254_RF1.mat',
% method of cwt, L parameter of -0.1254 and channel (option) 22
% then before saving, run the code below to add in other methods
clearvars -except dat channels thresholds
% fileName = '200127_FTDOrg_GrpD_5B_Slice7.mat' 
% fileName = '191210_slice1_DIV_g04_2018_2min.mat'
fileName = 'MPT190403_6B_DIV28.mat'
method = 'cwt'
parameter = -0.0627
refPeriod_ms = 1

if ~exist('dat','var')
    load(fileName)
end

elecpos = find(channels == 33)
trace = dat(:,elecpos);
% [spikeTrain, filtdata, threshold] = detectSpikes(trace,'cwt', 4.5,-0.188,refPeriod_ms);
[spikeTrain, filtdata, threshold] = detectSpikes(trace,'Manuel',3 ,-0.188,refPeriod_ms);
sum(spikeTrain)
sp_times=find(spikeTrain==1);
%%
n_spikes_to_plot=50;
%added correction if num spikes is fewer than desired number:
if  sum(spikeTrain) < n_spikes_to_plot
    n_spikes_to_plot = sum(spikeTrain);
end

time_win        = 30; % in ms

if ~exist('fs') %correction if fs isn't loaded properly (sampling freq)
    fs = 25000;
end

sample_win      = time_win * (fs/1000); % because sampling freq, fs, is in seconds
slidage         = 5 * fs/1000; % fs = 1 s; so this slides every 5 ms
num_windows     = (length(spikeTrain)- sample_win) / slidage;
num_done        = 0;
%progressbar
for i = 1:num_windows %-1 so that i can start from 0
    spike_sums(i)   = sum(spikeTrain(1+slidage*num_done:slidage*num_done+sample_win));
    num_done        = num_done + 1;
    %progressbar(i/num_windows)
end
most_spikes_index = find(spike_sums == max(spike_sums))-1; % -1 is necessary
%because if the first window had the most spikes, we need to start from
%one; if the 2nd window had the most spikes, we need to start from
%1+slidage*1

if  length(most_spikes_index) > 1
    most_spikes_index = most_spikes_index(2); %if it's a draw just take second one
else
end

%get index of samples to plot
samples_index = 1+slidage*most_spikes_index:slidage*most_spikes_index+sample_win;

%% plot
f1 = figure
x1 = (1:length(samples_index)) / (fs/1000); %time vec in ms
y1 = filtdata(samples_index);
plot(x1,y1,'Color','k');
aesthetics
box off
%ylim([min(y1) max(y1)])
ylim([-50 50])
hold on
%y2 = mean(y1)+6*std(y1);
y2 = 20;
sp_times = find(spikeTrain(samples_index) == 1);
% plot(x1(sp_times),y2,'v','MarkerSize',5,'MarkerFaceColor','b','MarkerEdgeColor','b');

removeAxis
sb=scalebar;
sb.Position=[0.5,-30];
sb.XLen = 2.5; % in ms
sb.YLen = 5; % in uV
sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
sb.hTextY_Pos= [-1000,-1000];

if      strcmp(method,'Manuel') | strcmp(method,'abs')
    threshvec = ones(1,length(x1))*threshold;
    hold on
%     plot(x1,threshvec,'LineStyle','--','Color','b','LineWidth',1.5)
    hold off
end

fileName1 = fileName;
if  contains(fileName,'_') %remove underscores for title
    fileName1(strfind(fileName1,'_'))=' ';
    fileName1=strcat('{',fileName1(1:end-4),'}');
    title({[int2str(length(samples_index)/fs*1000),'{ ms of recording from:}'],[fileName1],...
        ['{channel (MCS xy coordinate): }',num2str(channels(elecpos))],['{ (scalebars: horizontal }',num2str(sb.XLen),'{ ms}','{ vertical }',...
        int2str(sb.YLen),'{ uV}',')'],[' ']});
else
    title({[int2str(length(samples_index)/fs*1000),'{ ms of recording from:}'],[fileName1],...
        ['{channel (MCS xy coordinate): }',num2str(channels(electrode_to_plot(ei)))],['{ (scalebars: horizontal }',num2str(sb.XLen),'{ ms}','{ vertical }',...
        int2str(sb.YLen),'{ uV}',')'],[' ']});
end

%make title font smaller
ax = gca;
ax.TitleFontSizeMultiplier = 0.9;

%                 %% save as PNG
%                 fprintf(strcat('\n','\n',fileName(1:end-4),'\n','\n',' saving fig...', '\n','\n'))
%                 saveas(gcf,strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
%                     '_',num2str(time_win),'ms_trace_spikes_marked.png'));
%                 close all

% below corrections for organoid comparison
% ylim([-100 100])
% sb.Position=[0.5,-80];
% y2 = 30
%%
% [spikeTrain1, ~, ~] = detectSpikes(trace,'Manuel', 3,-0.4387,refPeriod_ms);
% sp_times1 = find(spikeTrain1(samples_index) == 1);
% p1 = plot(x1(sp_times),y2,'v','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 1],'MarkerEdgeColor',[0.5 0.5 1]);
p1 = plot(x1(sp_times),y2+5,'v','MarkerSize',5,'MarkerFaceColor',[0.8 0.6 0.6],'MarkerEdgeColor',[0.8 0.6 0.6]);
threshvec = ones(1,length(x1))*threshold;
plot(x1,threshvec,'LineStyle','--','Color',[0.8 0.6 0.6],'LineWidth',1.5)

% mecp2:
[spikeTrain2, ~, ~] = detectSpikes(trace,method, 0,-0.3761,refPeriod_ms);
% organoid:
% [spikeTrain2, ~, ~] = detectSpikes(trace,method, 0,-0.1880,refPeriod_ms);
sp_times2 = find(spikeTrain2(samples_index) == 1);
p2 = plot(x1(sp_times2),y2+15,'v','MarkerSize',5,'MarkerFaceColor',[0 0 0.5],'MarkerEdgeColor',[0 0 0.5]); beep

% mecp2:
[spikeTrain3, ~, threshold] = detectSpikes(trace,'cwt', 3,-0.5014,refPeriod_ms);beep
% organoid:
% [spikeTrain3, ~, threshold] = detectSpikes(trace,'cwt', 3,-0.2507,refPeriod_ms);beep
sp_times3 = find(spikeTrain3(samples_index) == 1);
p3 = plot(x1(sp_times3),y2+20,'v','MarkerSize',5,'MarkerFaceColor',[0 0 0.9],'MarkerEdgeColor',[0 0 0.9]);
% threshvec = ones(1,length(x1))*threshold;
% plot(x1,threshvec,'LineStyle','--','Color',[0.5 0 0],'LineWidth',1.5)

[spikeTrain4, ~, threshold] = detectSpikes(trace,'Manuel', 4.5,-0.3761,refPeriod_ms);
sp_times4 = find(spikeTrain4(samples_index) == 1);
try
    p4 = plot(x1(sp_times4),y2+10,'v','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
catch
end
threshvec = ones(1,length(x1))*threshold;
plot(x1,threshvec,'LineStyle','--','Color',[1 0 0],'LineWidth',1.5)

[spikeTrain5, ~, threshold] = detectSpikes(trace,'Manuel', 2.5,-0.3761,refPeriod_ms);
sp_times5 = find(spikeTrain5(samples_index) == 1);
p5 = plot(x1(sp_times5),y2,'v','MarkerSize',5,'MarkerFaceColor',[0.5 0 0],'MarkerEdgeColor',[0.5 0 0]);
threshvec = ones(1,length(x1))*threshold;
plot(x1,threshvec,'LineStyle','--','Color',[0.5 0 0],'LineWidth',1.5)


%%
saveas(gcf,strcat(fileName(1:end-4),'_channelXY_',int2str(channels(elecpos)),...
    '_',num2str(time_win),'ms_trace_compared_to_other_methods3.png'));




