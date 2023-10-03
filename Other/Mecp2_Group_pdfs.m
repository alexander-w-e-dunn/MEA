% Mecp2 figures
% code_dir = 'D:\MECP2_2019_AD\Scripts_and_Output';
code_dir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
addpath(genpath(code_dir));
cd(code_dir)
spikeDir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\Mecp2\Spikes';

%% plot adjM for cultures across genotypes at each DIV
% XLSdir = 'D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\Unsorted_byFileType\XLS_unsorted';
% PNGdir = 'D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\By_culture\MPT';
XLSdir          = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\XLS_unsorted';
AdMdir          = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\Mecp2\AdjMs';
% file_suffix     = '_DIV28_cSpikes_L-0.0627_RF1_adjM_0.05.mat';

cd(XLSdir)
genotypes   = readtable('files_included_gt.csv');
DIVs        = [14 21 28 35];

genotype_strings = flip(unique(genotypes.genotype))

cd(AdMdir)
onlyWT = genotypes(strcmp(genotypes.genotype,'WT'),:);
onlyKO = genotypes(strcmp(genotypes.genotype,'KO'),:);
onlyHE = genotypes(strcmp(genotypes.genotype,'HE'),:);
%% plot div14-div35 for each culture for a given genotype

for genotype = 1 : length(genotype_strings)
    f = figure;
    f.Position = [20.3333 41.6667 220.6667 599.3333];
    % get index of cultures that are this genotype
    %     for i = 1 : length(genotypes.culture)
    %        if strcmp(genotypes.genotype,genotype_strings(genotype))
    %     end
    num_culture = sum(strcmp(genotypes.genotype,genotype_strings(genotype)));
    cultureIDs  = find(strcmp(genotypes.genotype,genotype_strings(genotype)) == 1);
    row_num     = 0;
    plot_num    = 0;
    
    for culture = cultureIDs'
        row_num         = row_num+1;
        fileName        = genotypes{culture,1}{1};
        
        for age = DIVs
            file_suffix     = strcat('_DIV' ,num2str(age), '_cSpikes_L-0.0627_RF1_adjM_0.05.mat');
            plot_num = plot_num + 1;
            load(strcat(fileName,file_suffix));
            subplot(length(cultureIDs),length(DIVs),plot_num)
            imagesc(abs(adjM)); caxis([0 1]);
            
            % store adjM
            adjM_all{genotype}{row_num,find(age==DIVs)} = adjM;
            
            % remove y labels for all matrices not on left
            if age > 14
                yticklabels([])
            else %can remove this option
                % plot electrode label and arrow instead of IDs
                yticks(round(median(1:length(adjM))));
                yticklabels(['elec→']); ytickangle(90);
            end
            % remove x labels for all but last four plots
            if row_num < num_culture
                xticklabels([])
            else
                xticks(round(median(1:length(adjM))));
                xticklabels(['elec→']); %ytickangle(90);
            end
            % add div labels
            if row_num == 1
                title(strcat('DIV',num2str(age)))
            end
            
            clear adjM
            
        end
        
    end
    % save as pdf then close
    % saveas(f,strcat('all_adjMs_',genotype_strings{genotype},'.pdf'))
    f.PaperPositionMode='auto';
    print(strcat('all_adjMs_',genotype_strings{genotype},'.pdf'),...
        '-dpdf','-bestfit')
    close(f)
end

%% plot the mean adjM across cultures at each DIV for each genotype
%{
issue here is that edges don't match up between recordings so there is
dilution of edge weights across cultures. So this plot is not useful.
Instead get the vector of ABSOLUTE edge weights >0 across all recordings
and plot distribution for each GT over each DIV. Repeat for threshold of
0.1 (so distribution is not biased by inactive electrodes)
%}
f2 = figure;
plot_num = 0;
for genotype = 1 : length(genotype_strings)
    adjMs = adjM_all{genotype};
    for age = DIVs
        % adjMs_DIV = adjMs{:,find(age==DIVs)};
        for culture = 1 : size(adjMs,1)
            adjMs_DIV(:,:,culture) = adjMs{culture,find(age==DIVs)};
        end
        plot_num = plot_num + 1;
        subplot(length(genotype_strings) , length(DIVs) , plot_num )
        % change NaNs to count as a 0
        adjMs_DIV(find(isnan(adjMs_DIV)==1)) = 0;
        imagesc(nanmedian(abs(adjMs_DIV),3)); caxis([0 1]);
        % title(strcat(genotype_strings(genotype),'DIV',num2str(age)))
        if age > 14
            yticklabels([])
        else %can remove this option
            % plot electrode label and arrow instead of IDs
            yticks(round(median(1:length(adjM_all{1,1}{1,1}))));
            yticklabels(['elec→']); ytickangle(90);
            ylabel({[genotype_strings{genotype} '    ']});
            set(get(gca,'ylabel'),'rotation',0,'FontWeight','bold')
        end
        % remove x labels for all but last four plots
        if genotype < length(genotype_strings)
            xticklabels([])
        else
            xticks(round(median(1:length(adjM_all{1,1}{1,1}))));
            xticklabels(['elec→']); %ytickangle(90);
        end
        % add div labels
        if genotype == 1
            title(strcat('DIV',num2str(age)))
        end
    end
    clear adjMs adjMs_DIV
end

f2.PaperPositionMode='auto';
print(strcat('ave_adjMs_all_genotypes','.pdf'),...
    '-dpdf','-bestfit')
close(f2)

%% Raster plot over dev for each culture
% make page landscape and column for each DIV.

cd(spikeDir)

for genotype = 1 : length(genotype_strings)
    f = figure;
    f.Position = [20.3333 41.6667 886 599.3333];
    % get index of cultures that are this genotype
    %     for i = 1 : length(genotypes.culture)
    %        if strcmp(genotypes.genotype,genotype_strings(genotype))
    %     end
    num_culture = sum(strcmp(genotypes.genotype,genotype_strings(genotype)));
    cultureIDs  = find(strcmp(genotypes.genotype,genotype_strings(genotype)) == 1);
    row_num     = 0;
    plot_num    = 0;
    
    for culture = cultureIDs'
        row_num         = row_num+1;
        fileName        = genotypes{culture,1}{1};
        
        for age = DIVs
            file_suffix     = strcat('_DIV' ,num2str(age), '_cSpikes_L-0.0627_RF1.mat');
            plot_num = plot_num + 1;
            
            load(strcat(fileName,file_suffix));
            subplot(length(cultureIDs),length(DIVs),plot_num)
            % plot raster
            %remove spikes from ref
            spikeMatrix = zeros(size(cSpikes));
            spikeMatrix = full(cSpikes);
            spikeMatrix(:,find(channels==15))=zeros(size(spikeMatrix(:,find(channels==15))));
            
            if ~exist('fs')
                fs=25000;
            end
            
            if  floor(length(spikeMatrix)/fs)~=length(spikeMatrix)/fs;
                %calculate number of samples to subtract to make
                %length of rec in s a whole number
                n2del = fs*(length(spikeMatrix)/fs - floor(length(spikeMatrix)/fs));
                spikeMatrix=spikeMatrix(1:length(spikeMatrix)-(n2del),:);
            else
            end
            
            recordDuration = round(length(spikeMatrix)); %in samples
            %downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration * 1/fs);
            downFactor = 25000;
            downSpikeMatrix = downSampleSum(spikeMatrix, recordDuration/downFactor);
            new_fs = fs/downFactor;
            
            
            imagesc(log10(downSpikeMatrix'))
            %     imagesc(downSpikeMatrix')
            
            %aesthetics
            %ylabel('Electrode')
            %xlabel('Time (s)')
            %cb = colorbar;
            % ylabel(cb, 'Spike count')
            %ylabel(cb, 'Firing Rate (Hz)')
            %cb.TickDirection = 'out';
            % cb.Ticks = 0:5; % for slice 5 specifically
            %set(gca,'TickDir','out');
            %cb.Location = 'Southoutside';
            %cb.Box = 'off';
            %set(gca, 'FontSize', 14)
            ylimit_cbar = 3;
            caxis([0,ylimit_cbar]) %set colorbar axis limits; also adjusts colour
            %below command does not adjust colour hence need for caxis command
            %remove decimal ticks e.g. 0.5 Hz
            %cb.Ticks = linspace(0,ylimit_cbar,ylimit_cbar/1+1);%(start,end,number of numbers)
            %below is for plotting log scale, labelling in raw spike count
            %cb.TickLabels = 10.^(linspace(0,ylimit_cbar,ylimit_cbar+1));
            
            %yticks([1, 10:10:60])
            
            
            
            
            % imagesc(abs(adjM)); caxis([0 1]);
            
            % store adjM
            % adjM_all{genotype}{row_num,find(age==DIVs)} = adjM;
            
            % remove y labels for all matrices not on left
            if age > 14
                yticklabels([])
            else %can remove this option
                % plot electrode label and arrow instead of IDs
                yticks(round(median(1:size(downSpikeMatrix,2))));
                yticklabels(['Elec→']); ytickangle(90);
            end
            % remove x labels for all but last four plots
            if row_num < num_culture
                xticklabels([])
            else
                xticks([1 round(median(1:length(downSpikeMatrix))) length(downSpikeMatrix)]);
                xticklabels({'0','Time (s)→',num2str(length(downSpikeMatrix))}); %ytickangle(90);
            end
            % add div labels
            if row_num == 1
                title(strcat('DIV',num2str(age)))
            end
            
            clear downSpikeMatrix 
            spikeMatrix = zeros(size(cSpikes));
            
        end
        
    end
    % save as pdf then close
    % saveas(f,strcat('all_adjMs_',genotype_strings{genotype},'.pdf'))
    f.PaperOrientation = 'landscape'; f.PaperPositionMode='auto'; 
    print(strcat('all_Rasters_',genotype_strings{genotype},'.pdf'),...
        '-dpdf','-bestfit')
    close(f)
end