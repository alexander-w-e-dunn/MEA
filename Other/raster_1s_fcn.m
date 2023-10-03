% Plot a raster of an MEA
% currently only 1 s plot with ms axis or scale bar fully supported
% need to add support for seconds scale bar
% add inputs scalebarVSaxis; timescale of axes; length of rec. to plot
% 
% INPUTS
%   spikeMatrix:    num samples x num channels, binary matrix, 1 = spike
%   fs:             sampling frequency
% 
% to do:
% add the option to plot the spikiest 1 s of activity rather than the first
% second, this may be useful for comparing before and after TTX
% 
% Author: AWE Dunn; U. of Cambrirdge; June 2020
% Last edit: 3 June 2020

function raster_1s_fcn(spikeMatrix,fs)

%remove spikes from ref
try
    spikeMatrix(:,find(channels==15))=zeros(size(spikeMatrix(:,find(channels==15))));
catch
    disp('no channels variable cannot ground elec')
end

spikeMatBackup = spikeMatrix;
s_to_plot = 1;
if ~exist('fs','var')
    fs = 25000;
end
spikeMatrix=spikeMatrix(1:fs*s_to_plot,:);

%% downsample to 1000 hz

if fs == 25000
    
    while  floor(length(spikeMatrix)/fs)~=length(spikeMatrix)/fs;
        %calculate number of samples to subtract to make
        %length of rec in s a whole number
        n2del = fs*(length(spikeMatrix)/fs - floor(length(spikeMatrix)/fs));
        spikeMatrix=spikeMatrix(1:length(spikeMatrix)-(n2del),:);
    end
    
    recordDuration = round(length(spikeMatrix)); %in samples
    %downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration * 1/fs);
    downFactor = 25;
    if ~~exist('rand_spikeMat_sparse_all')
        downFactor = 1000;
    end
    
    spikeMatrix = downSampleSum(spikeMatrix, recordDuration/downFactor);
    fs = fs/downFactor;
    
end

spikeMatrix(find(spikeMatrix > 1)) = 1;

% figure
% clear T
%     T = find(spikeMatrix(:) == 1);
%     N = size(spikeMatrix,2);
%     L = length(spikeMatrix);
%     %   RASTERPLOT(T,N,L) Plots the rasters of spiketimes (T in samples) for N trials, each of length
%     %   L samples, Sampling rate = 1kHz. Spiketimes are hashed by the trial length.%
%     %   RASTERPLOT(T,N,L,H) Plots the rasters in the axis handle H%
%     %   RASTERPLOT(T,N,L,H,FS) Plots the rasters in the axis handle H. Uses sampling rate of FS (Hz)
%     rasterplot(T(elec),N,L) % (askrajiv@gmail.com; Boston University, Boston, MA)

%% plot
progressbar('electrodes')
lineshaving = 0.15;
f1 = figure; 
f1.Position = [680 50 560 800];
axisscalebar = 'scalebar';
timescale = 'milliseconds';
    
for elec = 1 : size(spikeMatrix,2)
    
    spiketimes = find(spikeMatrix(:,elec));
    
    ypos = elec - 1;
   
    ybot = zeros(length(spikeMatrix),1) + ypos + lineshaving;
    ytop = ybot;
    ytop(spiketimes) = ypos + 1 - lineshaving;
    
    xx = [1:length(spikeMatrix) ; 1:length(spikeMatrix)];
    yy = [ybot' ; ytop'];
    subplot(5,1,[2 3 4 5])
    p1 = plot(xx, yy, 'k'); 
    if elec == 1
        ylim([0 size(spikeMatrix,2)-1])
        
        aesthetics
        hold on
        set(gca,'TickDir','out','LineWidth',2)
        ylabel('electrode');
        xlabel('time (ms)');
        
        if strcmp(timescale,'seconds')   ;
            xlabel('time (s)');
            xticklabels(xticks ./ fs);
        else
            xticklabels(xticks ./ (fs/1000) );
        end
        
        set(gca,'FontName','Arial','FontSize',14)
        
        if strcmp(axisscalebar,'scalebar')
            % currently only ms is supported
            axis off
            sb = scalebar;
            ylim([-1 size(spikeMatrix,2)]);
            sb.Position = [1 -0.5];
            sb.YLen = 0;
            sb.hTextY.String = '';
            sb.XLen = 100;
            sb.hTextX.FontName = 'Arial';
            sb.hTextX.FontSize = 14;
            sb.XUnit = 'ms';
        end
        
    end
    clear spiketime ypos ytop ybot xx yy
    %     disp('elec done')
    progressbar(elec/size(spikeMatrix,2))
end

%% add histogram of 10 ms time bins of spikes across whole array
time_bin_ms = 10;
spikecounts_10msBins = sum(downSampleSum(spikeMatrix,length(spikeMatrix)/time_bin_ms),2);
x = linspace(1,length(spikeMatrix),length(spikecounts_10msBins));

subplot(5,1,1)
bar(x,spikecounts_10msBins,'k') 
aesthetics
axis off
ylim([0 50])
sb2 = scalebar;
sb2.Position = [1 30];
sb2.XLen = 0; sb2.hTextX.String = '';
sb2.YLen = 25;
sb2.hTextY.FontName = 'Arial';
sb2.hTextY.FontSize = 14;
sb2.hTextY.Position = [-45  31.0000         0];
sb2.YUnit = 'spikes';

end
