function batch_getRaster_fcn(files)

%code:
progressbar('files completed')
for i=1:length(files)
    progressbar(i/length(files));
    fileName    =   files(i).name;
    if ~exist(strcat(fileName(1:end-4),'_Raster.png'))
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
            
        elseif ~~exist('rand_spikeMat_sparse_all')
            spikeMatrix = full(rand_spikeMat_sparse_all{1,1});
            fs = 1000;
        else
            disp('error - cant get spike mat')
        end
        
        %remove spikes from ref
        try
            spikeMatrix(:,find(channels==15))=zeros(size(spikeMatrix(:,find(channels==15))));
        catch
            disp('no channels variable cannot ground elec')
        end
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
        
        while  floor(length(spikeMatrix)/fs)~=length(spikeMatrix)/fs;
            %calculate number of samples to subtract to make
            %length of rec in s a whole number
            n2del = fs*(length(spikeMatrix)/fs - floor(length(spikeMatrix)/fs));
            spikeMatrix=spikeMatrix(1:length(spikeMatrix)-(n2del),:);
        end
        
        recordDuration = round(length(spikeMatrix)); %in samples
        %downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration * 1/fs);
        downFactor = 25000;
        if ~~exist('rand_spikeMat_sparse_all')
            downFactor = 1000;
        end
        
        downSpikeMatrix = downSampleSum(spikeMatrix, recordDuration/downFactor);
        new_fs = fs/downFactor;
        
        % reorder by columns, ascending order
        try 
            [~, I] = fcn_sort_channels_MEA(channels,'columns');
            downSpikeMatrix = downSpikeMatrix(:,I);
        catch
        end
        
        figure1 = figure;
        %imagesc(log10(downSpikeMatrix'))
        downSpikeMatrix(:,find(sum(downSpikeMatrix)==0)) = NaN;
        h = imagesc(log10(downSpikeMatrix'));
        
        % black background - currently doesnt work
        %         axes1 = axes('Parent', figure1); hold(axes1,'on');
        %         rectangle('Position',[0.5000 0.5000 720.5000 60.5000],'FaceColor',[0 0 0]);
        %         hold on; h = imagesc(downSpikeMatrix'); hold off;
        
        %         % make the quiescent periods white & make removed channels white
        downSpikeMatrix(find(downSpikeMatrix==0)) = NaN;
        set(h, 'AlphaData', ~isnan(downSpikeMatrix'));
        
        % make the quiescent periods and removed channels black
        %           c = parula; c(1,:) = [0 0 0]; colormap(c); %also adds black section as 0 value of colorbar
        
        % change color scheme
        %         c = bone; c = c(1:200,:); c = flipud(c);
        %         c = copper;
        %         c = jet; c = c(end - round(length(c)/3) : end , :); %c = flipud(c);
        %         c = parula; c = c(end - round(length(c)/2) : end , :); c = flipud(c);
                c = parula;c = c(1:round(length(c)*.85),:);
%         c = gray; c = c(1:250,:); c = flipud(c);
        colormap(c);
        
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
%         caxis([0,ylimit_cbar]) %set colorbar axis limits; also adjusts colour
        caxis([-1,ylimit_cbar])
        %below command does not adjust colour hence need for caxis command
        %remove decimal ticks e.g. 0.5 Hz
%         cb.Ticks = linspace(0,ylimit_cbar,ylimit_cbar/1+1);%(start,end,number of numbers)
        cb.Ticks = log10([0.1 0.2 0.5 1 2 5 10 20 50 100]);
        %below is for plotting log scale, labelling in raw spike count
%         cb.TickLabels = 10.^(linspace(0,ylimit_cbar,ylimit_cbar+1));
        cb.TickLabels = 10 .^ cb.Ticks;

        yticks([1, 10:10:60])
        
        fileName1=files(i).name;
        if  contains(files(i).name,'_') %remove underscores for title
            fileName1(strfind(fileName1,'_'))=' ';
            %fileName1=strcat('{',fileName1,'}');
            title({strcat(fileName1(1:end-4),' Raster'),' '});
        else
            title({strcat(files(i).name(1:end-4),' Raster'),' '});
        end
        %make title font smaller
        ax = gca;
        ax.TitleFontSizeMultiplier = 0.7;
        
        %save raster as PNG
        fprintf(strcat('\n','\n',files(i).name(1:end-4),' saving raster...', '\n','\n'))
        saveas(gcf,strcat(fileName(1:end-4),'_Raster.png'));
        if ~~exist('figure1')
            close(figure1); 
        end
        clear cSpikes mSpikes aSpikes spikeMatrix fs downSpikeMatrix fileName fileName1
    else
        disp('file done')
    end
end


end
