% get p dist of getting 20 spikes in any 1 s window across range of FRs
%  get(groot,'default')

set(groot,'defaultAxesFontName','Arial')
set(groot,'defaultAxesFontSize',22)
set(groot,'defaultAxesLineWidth',2)
set(groot,'defaultFigurePosition',[700   435   520   350])

f1 = figure;
lambda = 23;
FR = [0:0.01:40];
r = poissrnd(lambda,[10000 1]);
numberOfBins = length(unique(r));
h = histfit(r,numberOfBins,'normal'); aesthetics
h(1).FaceColor = 0.5*[1 1 1]; h(1).EdgeColor = 1*[1 1 1]; h(1).LineWidth = 2;
h(2).LineWidth = 2; h(2).Color = 0*[1 1 1]; h(2).LineStyle = '--';
yt = get(gca, 'YTick'); yt2 = yt/numel(r);
ytlabels = yt/numel(r);
set(gca, 'YTick', yt, 'YTickLabel', ytlabels);
xlabel('spike count'); ylabel('p');
% l1 = legend('Box','off'); 
% l1.Title.String = {'\mu{\itr} (Hz)'}; 
% l1.Title.FontWeight = 'normal';
realxlims = xlim; realylims = ylim;

saveas(f1,'p of n spikes in 1 s time bins.png')



f2 = figure;
count = 0; linestyles = {'-', '--', '-.',':'};
for lambda = [1 5 10 20]
    pd = makedist('Poisson','lambda',lambda);
    y = pdf(pd,1:40);
    count = count+1;
    plot(y,'k','LineStyle',linestyles{count}); aesthetics;
    xlabel('spike count'); ylabel('p');
    hold on
end
l2 = legend('Box','off');
l2.String = {'1' '5' '10' '20'};
l2.Title.String = '\mu{\itr} (Hz)';
l2.Title.FontWeight = 'normal';
ylims = ylim;
ylim([ylims(1)-0.01 ylims(2)])
saveas(f2,'p of n spikes in 1 s time bins for each rate pdf.png')


fileName = '200127_FTDOrg_GrpD_5B_Slice7_mSpikes_3.mat'
load(fileName)
spikeMatrix = full(mSpikes);
spikecounts = downSampleSum(spikeMatrix(:,1),360);
f3 = figure;
numberOfBins = length(unique(spikecounts));
h = histfit(spikecounts,numberOfBins,'normal'); aesthetics
h(1).FaceColor = 0.5*[1 1 1]; h(1).EdgeColor = 1*[1 1 1]; h(1).LineWidth = 2;
h(2).LineWidth = 2; h(2).Color = 0*[1 1 1]; h(2).LineStyle = '--';
yticks(yt2 * numel(spikecounts)); ylim([min(yticks) max(yticks)]);
set(gca, 'YTick', yt, 'YTickLabel', ytlabels);
xlabel('spike count'); ylabel('p');
xlim(realxlims); 
yticks(yt2 * numel(spikecounts));
saveas(f3,'actual proportion of n spikes in 1 s time bins.png')

