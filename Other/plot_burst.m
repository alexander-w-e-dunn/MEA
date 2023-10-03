% plot burst
clear all
load('burst2plot.mat')
load(fileName)
% get only active channels
dat = dat(burst_time(1):burst_time(2),burstchans);
% filter and downsample
for elec = 1:length(burstchans)
    [~, filtdata(:,elec), ~] = detectSpikes(dat(:,elec),'Manuel', 3,-0.188,1.5);
end
t2plot_ms = 25;
kHz = 10;
fsms = fs/1000;
new_fs = kHz * 1000;
new_n_samples = t2plot_ms * kHz;
while  round(floor(length(filtdata)/fsms) , -1 ) ~= length(filtdata)/fsms;
    n2del = fsms * (length(filtdata)/fsms - round(floor(length(filtdata)/fsms) , -1 ));
    filtdata=filtdata(1:length(filtdata)-(n2del),:);
end
downfiltdat = downSampleMean(filtdata, new_n_samples);
MEAgraphics(1)
% plot(downfiltdat(:,1))
f1 = figure;
f1.Position = [ 228         143        1065         805];
% subplot(length(burstchans),1,1)
numchans = 6;
subplot(numchans+1,1,1)
c=lines;
for elec = 1:numchans %length(burstchans)
%     subplot(length(burstchans),1,elec)
    subplot(numchans+1,1,elec)
    plot(downfiltdat(:,elec),'LineWidth',2,'Color', c(elec,:) )
    axis off
end
% subplot(numchans+1,1,4)
% linkaxes
axis off
% plot spike raster underneath
% could plot 10 ms before the onset of the burst
% then indicate onset of burst with a line?
% chans = [4 5 6];
% stimes1 = find( burst( 1:t2plot_ms*fsms , burstchans(4) ) == 1 );
for elec = 1 : numchans
    m = 1.75;
    thr = median(downfiltdat(:,elec)) - (m * std(downfiltdat(:,elec))) ;
    stimes{:,elec} = find(downfiltdat(:,elec) < thr) 
end 


for elec = 1 : numchans
    
    subplot(numchans+1,1,numchans+1)
    plot([stimes{:,elec} , stimes{:,elec}],[0  1],...
        'Color',c(elec,:) ) ; 
    hold on
end 
aesthetics
box off; axis off

subplot(numchans+1,1,numchans+1)
sb = scalebar;
sb.Position = [10 -1];
sb.YLen = 0.5;
sb.hTextX_Pos = [-10 -10]; 
sb.hTextY_Pos = [-10 -10]; 


