close all

for spike = 1 : 50
    plot(spikeWaveforms{1,1}.mea(spike,:))
    hold on
end

plot(mean(spikeWaveforms{1,1}.mea(:,:)),'-r','LineWidth',2)

figure
for spike = 1 : 50
    plot(spikeWaveforms{1,1}.thr3p5(spike,:))
    hold on
end

plot(mean(spikeWaveforms{1,1}.thr3p5(:,:)),'-r','LineWidth',2)


%% all methods
scriptsDir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
addpath(genpath(scriptsDir));
kHz = 25;

methods = {'mea';'bior1p5';'thr3p5';'thr5'};

for method = 1 : length (methods)
    for spike = 1 : 50
        meanWaves{method} = mean(spikeWaveforms{1,10}.(methods{method})(:,:));
    end
end

figure
MEAgraphics(1)
for method = 1 : length (methods)
    p = plot(meanWaves{method},'LineWidth',1.5)
    p.Color(4) = 0.5
    hold on
end
legend(methods)
aesthetics
xticklabels(xticks / kHz);
xlabel('Time (s)'); ylabel('\muV');

%% compare mea and bior1p5
methods = {'mea';'bior1p5'};

for method = 1: length(methods)
    for channel = 1:length(spikeTimes)
        spikeCounts(method,channel) = length(spikeTimes{channel}.(methods{method})(:,:));
    end
end

figure
scatter(spikeCounts(1,:),spikeCounts(2,:),'LineWidth',2)
title('Spike counts for each electrode per method')
xlabel('mea'); ylabel('bior1.5');
hold on; x = xlim;
plot(1:x(2),':k')
scatter(spikeCounts(1,:),spikeCounts(2,:),'LineWidth',2)
