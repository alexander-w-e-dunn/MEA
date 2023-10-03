% MAKE 100 s raster
refPeriod_ms = 1;
close all; clearvars -except refPeriod_ms
fprintf(strcat('\n','\n','creating raster plots...','\n','\n'))

%inputs:
% files = dir('*FTD*mSpikes_3.mat');                   % where your .mat files are
files = dir('MPT190403*6B*cSpikes_L0*.mat');                   % where your .mat files are
files = files(~~contains({files.name}, 'TTX'));
% files = files(~contains({files.name}, 'ttx'));
% files = files(~contains({files.name}, 'stim'));
% files = files(~contains({files.name}, 'edited'));
% files = files(~contains({files.name}, '2min'));
% files = files(~contains({files.name}, 'adjM'));
% files = files(~contains({files.name}, 'DIV07'));
% files = files(~~contains({files.name}, '200107') & contains({files.name}, 'DIV21'));
% %files = files(~contains({files.name}, '_4'));       % remove ones already done
% % files = files(41:end);

%code:
for i=1:length(files)
    fileName    =   files(i).name
    load(fileName)
    if ~exist('fs')
        fs=25000;
    else
    end
    
    if ~~exist('aSpikes')
        spikeMatrix = full(aSpikes);
        
    elseif ~~exist('cSpikes')
        spikeMatrix = full(cSpikes);
        
    elseif ~~exist('mSpikes')
        spikeMatrix = full(mSpikes);
        
    else
        disp('error - cant get spike mat')
    end
    
    %remove spikes from ref
    spikeMatrix(:,find(channels==15))=zeros(size(spikeMatrix(:,find(channels==15))));
    %sum(spikeMatrix(:,find(channels==15))); %check no spikes
    %correction if length rec in s is not a whole number
    
    %    if ~isinteger(length(spikeMatrix)/fs)
    %        %calculate number of samples to subtract to make
    %        %length of rec in s a whole number
    %        n2del = fs*(length(spikeMatrix)/fs - round(length(spikeMatrix)/fs));
    %        spikeMatrix=spikeMatrix(1:length(spikeMatrix)-n2del,:);
    %    else
    %    end
    
    %recordDuration = length(spikeMatrix); %in samples
    %downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration * 1/fs);
    
    if  floor(length(spikeMatrix)/fs)~=length(spikeMatrix)/fs;
        %calculate number of samples to subtract to make
        %length of rec in s a whole number
        n2del = fs*(length(spikeMatrix)/fs - floor(length(spikeMatrix)/fs));
        spikeMatrix=spikeMatrix(1:length(spikeMatrix)-(n2del),:);
    else
    end
    
    recordDuration = round(length(spikeMatrix)); %in samples
    %downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration * 1/fs);
    bin_size_s = 100; %bin size in s
    downFactor = bin_size_s * 25000;
    
    % make new number of samples a whole number
    if  floor(recordDuration/downFactor) ~= recordDuration/downFactor
        n2del = fs*(recordDuration/downFactor * bin_size_s - floor(recordDuration/downFactor) * bin_size_s);
        spikeMatrix=spikeMatrix(1:length(spikeMatrix)-(n2del),:);
    else
    end
    
    downSpikeMatrix = downSampleSum(spikeMatrix, floor(recordDuration/downFactor));
    new_fs = fs/downFactor;
    
    figure
    imagesc(log10(downSpikeMatrix'*new_fs))
%     imagesc(downSpikeMatrix')
    
    aesthetics
    ylabel('Electrode')
    xlabel('Time (s)')
    cb = colorbar;
    % ylabel(cb, 'Spike count')
    ylabel(cb, 'Firing Rate (Hz)')
    cb.TickDirection = 'out';
    % cb.Ticks = 0:5; % for slice 5 specifically
    set(gca,'TickDir','out');
    cb.Location = 'Southoutside';
    cb.Box = 'off';
    set(gca, 'FontSize', 14)
    ylimit_cbar = 2;
    ymin_cbar   = -2;
    caxis([ymin_cbar,ylimit_cbar]) %set colorbar axis limits; also adjusts colour
    cb.Ticks = linspace(ymin_cbar,ylimit_cbar,(ylimit_cbar-ymin_cbar)+1);
    cb.TickLabels = 10.^(linspace(ymin_cbar,ylimit_cbar,(ylimit_cbar-ymin_cbar)+1));
    
    if contains(fileName,'org','ignorecase',true)
        ylimit_cbar = 3;
        caxis([0,ylimit_cbar])
        cb.Ticks = linspace(0,ylimit_cbar,ylimit_cbar/1+1);
        cb.TickLabels = 10.^(linspace(0,ylimit_cbar,ylimit_cbar+1));
    end

    
    yticks([1, 10:10:60])
    
    xticklabels(xticks*bin_size_s);
    
    fileName1=files(i).name;
    if  contains(files(i).name,'_') %remove underscores for title
        fileName1(strfind(fileName1,'_'))=' ';
        %fileName1=strcat('{',fileName1,'}');
        title({strcat(fileName1(1:end-4),' ',num2str(bin_size_s),' s Raster'),' '});
    else
        title({strcat(files(i).name(1:end-4),' ',num2str(bin_size_s),' s Raster'),' '});
    end
    %make title font smaller
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.7;
    
    %save raster as PNG
    fprintf(strcat('\n','\n',files(i).name(1:end-4),' saving raster...', '\n','\n'))
    saveas(gcf,strcat(fileName(1:end-4),'_',num2str(bin_size_s),'_s_Raster.png'));
    close all; clear cSpikes mSpikes aSpikes spikeMatrix fs downSpikeMatrix fileName fileName1
end