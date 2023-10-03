%  script for FENS 2020 poster
% plot raw then filtered trace with scale bar and linked axes
clear all
prep_plots % set defaults for plotting 

% fileName = '200127_FTDOrg_GrpD_5B_Slice7.mat'
fileName = 'MPT190403_6B_DIV50TTXBEFORE.mat'
t2plot_ms = 50;

if ~exist('data','var')
    data = load(fileName)
    fs = data.fs;
end

kHz = 5; % downsample to 10 kHz
new_n_samples = t2plot_ms * kHz;
rawdat = data.dat(50+[1:t2plot_ms*25],19);
downrawdat = downSampleMean(rawdat, new_n_samples);
f1 = figure
subplot(2,1,1)
plot(downrawdat - mean(downrawdat));aesthetics;axis off

[spikes, filtdat, ~] = detectSpikes(rawdat,'abs',-20,0,1.5);
[spikes2, ~, ~] = detectSpikes(rawdat,'cwt',1,0,1);
length(find(spikes==1))
length(find(spikes2==1))
downfiltdat = downSampleMean(filtdat, new_n_samples);

subplot(2,1,2)
plot(downfiltdat);aesthetics; linkaxes; axis off
f1.Position = [700   435   260   350];
xlim([0 150])
% subplot(9,1,9)
% plot(downfiltdat);aesthetics; linkaxes; axis off
sb = scalebar;
sb.XLen = 2 * length(downfiltdat) / t2plot_ms; % 2 ms scale bar
sb.YLen = 15; % 15 uV scale bar
sb.hTextX_Pos = 1000 * sb.hTextX_Pos;sb.hTextY_Pos = 1000 * sb.hTextY_Pos;
sb.Position = [1 -50];
saveas(f1,'filtering.png')

%% apply spike det
f2 = figure;
subplot(2,3,1:3)
plot(downfiltdat);aesthetics; linkaxes; axis off
f2.Position = [700   435   260   350];
xlim([0 150])
% subplot(9,1,9)
% plot(downfiltdat);aesthetics; linkaxes; axis off
sb = scalebar;
sb.XLen = 2 * length(downfiltdat) / t2plot_ms; % 2 ms scale bar
sb.YLen = 15; % 15 uV scale bar
sb.hTextX_Pos = 1000 * sb.hTextX_Pos;sb.hTextY_Pos = 1000 * sb.hTextY_Pos;
sb.Position = [1 -50];
xvals = xlim;
x = xvals(1):xvals(2);
plot(x,-20*ones(size(x)),'r--');
threshspikes = find(spikes == 1) / kHz;
yvals = ylim; hold on
plot([threshspikes threshspikes]', (ones(size(threshspikes)) * [yvals(2)-8 yvals(2)])','k-');
% saveas(f1,'filtering.png')
% figure
chans2plot = [19 16 24];
for i = 1:3
    subplot(2,3,i+3)
    [spikeTrain, finalData, ~] = detectSpikes(data.dat(:,chans2plot(i)),'abs',-20,0,1.5);
    spike_overlaysingle_fcn(spikeTrain,finalData,'abs',-20,1.5);
    if i == 1
        sb2=scalebar;
        sb2.XLen = 25;sb2.YLen = 15;
        sb2.hTextX_Pos = 1000 * sb2.hTextX_Pos;sb2.hTextY_Pos = 1000 * sb2.hTextY_Pos;
        sb2.Position = [-25 -40];
    end
end
saveas(f2,'applyingInitialParams.png')
close(f2)








%% 
% f2 = figure;
% h = histfit(filtdat,50,'normal');
% h(2).LineWidth = 3
% ranks = sort(h(2).XData);
% ranks(1);
% 
% alldat = data.dat(:,1);
% [~, allfiltdat, ~] = detectSpikes(alldat,'Manuel',5,0,1);
% t2plot_ms = length(allfiltdat)/fs*1000;
% % kHz = 5; % downsample to 10 kHz
% % new_n_samples = t2plot_ms * kHz;
% % downallfiltdat = downSampleMean(allfiltdat, new_n_samples);
% 
% f2 = figure;
% h = histfit(allfiltdat,1000,'normal');
% h(2).LineWidth = 3
% yt = get(gca, 'YTick'); 
% ytlabels = yt/numel(allfiltdat);
% set(gca, 'YTick', yt, 'YTickLabel', round(ytlabels,2));
% xlabel('\muV'); ylabel('p');aesthetics
% % h(1).FaceColor = 0.5*[1 1 1]; h(1).EdgeColor = 1*[1 1 1]; h(1).LineWidth = 2;
% % h(2).LineWidth = 2; h(2).Color = 0*[1 1 1]; h(2).LineStyle = '--';
% ranks = sort(allfiltdat);
% 
% 
% f3=figure;
% % get top 0.01 % of positive amps
% posranks = ranks(ranks > 0); posranks = posranks(posranks > prctile(posranks,99.99));
% % get that many neg amps
% negranks = ranks(1:length(posranks));
% % negranks = ranks(ranks < 0); negranks = negranks(negranks < prctile(negranks,99.99));
% % h1 = histfit(abs(ranks([1:500])),50,'normal')  
% % hold on; h2 = histfit(ranks([end-499: end]),50,'normal') 
% h1 = histfit(abs(negranks),50,'ev')  ;
% hold on; h2 = histfit(posranks,50,'ev') ;
% h1(2).Color = 'b';
% % h1(1).Visible = 'off'; 
% h2(1).Visible = 'off';
% % yt = get(gca, 'YTick'); 
% % ytlabels = yt/numel(allfiltdat);
% % set(gca, 'YTick', yt, 'YTickLabel', round(ytlabels,7));
% 
% 
% f4=figure;
% h1 = histfit(ranks,1000000,'normal')  ;
% h1(2).Color = 'b';
% % yt = get(gca, 'YTick'); 
% % ytlabels = yt/numel(allfiltdat);
% % set(gca, 'YTick', yt, 'YTickLabel', round(ytlabels,6));
% aesthetics
% ylabel('frequency'); xlabel('\muV')
% h1(1).FaceColor = 0.5*[0 0 1];% h1(1).EdgeColor = 0.5*[1 1 1]; h1(1).LineWidth = 1;
% h1(2).LineWidth = 2; h1(2).Color = 0*[1 1 1]; h1(2).LineStyle = '--';
% saveas(f4,'hist of amps all.png')
% xlims = xlim;
% xlim([xlims(1) -13])
% xvals = xticks;
% saveas(f4,'hist of amps neg.png')
% xlim([13 abs(xlims(1))])
% xticks(sort(abs(xvals)))
% saveas(f4,'hist of amps pos.png')


