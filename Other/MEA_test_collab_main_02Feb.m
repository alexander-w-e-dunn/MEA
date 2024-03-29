% analysis for MEA testing collaboration with Andras, Sagnik, Ilaria

%{

Note: exclude files containing stim artefacts and first run
remove_stim_artefacts_auto.m to get the cleaned files. Then you can run
this script on the cleaned files

%}

clear all
%% get files and save spikes

% inputs:

%{
spike detection methods: 'Manuel', 'cwt' or 'abs'
recommended parameter: 5 or 6, 0 and -20 to -25, respectively

note if spike detection method is 'abs' (absolute threshold), set the
'multiplier variable to be the threshold that you want (e.g. -20) and use
the minus sign if using a negative threshold
%}

% cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.1.File_Conversion_Scripts'
% MEAbatchConvert_alex
% cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'

% DATA AND FILES
clear all; close all
my_email    = 'awed2@cam.ac.uk';
scriptDir   = 'C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
dataDir     = 'F:\Alex\MAT files';
addpath(genpath(scriptDir))
cd(dataDir)
% data_and_scripts_dir = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
% cd(data_and_scripts_dir)
files = dir('MPA*.mat');                   
files = files(~contains({files.name}, 'TTX','IgnoreCase',true));
files = files(~contains({files.name}, 'Spikes'));
files = files(~contains({files.name}, 'stim'));

%SPIKE DETECTION METHODS AND PARAMETERS
%option one: two methods, three parameters for each
% meths   ={'Manuel';'abs'};
% params  =[3.5,3.25,3.75;                                % Threshold
%     -17,-18,-19]';                         % L parameter

%option two: two methods, one parameter for each
%meths   ={'Manuel';'cwt'};
%params  =[5;0]';                               % threshold; L parameter

%option three: one method, one parameter
meths   ={'cwt'};                              % one method only for speed
params  =[3,-0.0627];

%option four: one method, three parameters
% meths   ={'Manuel'};
% params  =[4,5,6;                                % Threshold
%     -0.1254,0,0.1254]';                         % L parameter

% custom
% meths   ={'Manuel';'cwt'};
% meths   ={'cwt'};
% params(:,1) = 3.2:0.2:5;
% params(:,2) = -log([10^0 10^1 10^2 10^3 10^4 10^5 10^6 10^7 10^8 10^9])/36.7368;
% params(:,1) = [2.5 3 3.5 4 4.5 5];
% params(:,2) = -log([10^0 10^1 10^2 10^3 10^4 10^6])/36.7368;

% CODE TO GET SPIKES AND SAVE:
fprintf(strcat('\n','\n','getting spike matrices...','\n','\n'))
disp(strcat({'number of files to do = '},num2str(length(files))))
progressbar('files','parameters','methods')

warning('off','MATLAB:load:variableNotFound') %turn warning off for now -
%just warns that a variable that i dont use anymore is not found

refPeriod_ms = 0; %choose refractory period in ms
for m = 1:length(meths)
    method              =       meths{m}
    
    for p = 1:size(params,1);
        L               =       params(p+size(params,1))
        multiplier      =       params(p)
        if strcmp(meths{m},'abs')
            multiplier  =       params(p+size(params,1))
            disp('adjusted multiplier for abs method')
        end
        
        batchGetSpike_function(dataDir,files,method,multiplier,L,refPeriod_ms,my_email)
        progressbar([],p/size(params,1),[]);       %update parameter progress
    end
    progressbar([],[],m/length(meths));            %update % methods done
end
warning('on','MATLAB:load:variableNotFound')

% outputStats_organoid
% disp(' comment out line 25�27 and 89�90 from MEA test collab script!!')
%%%%%%%%%%%%%%%%%%%%%%%%%
%% detect spikes in ttx file using corresponding none-ttx file thresholds in each channel
%%%%%%%%%%%%%%%%%%%%%%%%%
% for m = 1:length(meths)
%     method              =       meths{m}
%
%     for p = 1:size(params,1);
%         L               =       params(p+size(params,1))
%         multiplier      =       params(p)
%
%         batchGetSpike_function_TTX(data_and_scripts_dir,files,method,multiplier,L)
%         progressbar([],p/size(params,1),[]);       %update parameter progress
%     end
%     progressbar([],[],m/length(meths));            %update % methods done
% end

%% creater rasters and save
close all; clearvars -except refPeriod_ms
fprintf(strcat('\n','\n','creating raster plots...','\n','\n'))

%inputs:
files = dir('*FTD*slice*mSpikes_3.mat');  % where your .mat files are
% files = dir('*FTD*mSpikes_3.mat');                   % where your .mat files are
% files = dir('*MPT2001*cSpikes_L0.mat');                   % where your .mat files are
% files = files(~contains({files.name}, 'TTX'));
% files = files(~contains({files.name}, 'ttx'));
files = files(~~contains({files.name}, 'FTD'));
files = files(~contains({files.name}, '190814'));
files = files(17);
% % files = files(~contains({files.name}, '2min'));
% % files = files(~contains({files.name}, 'adjM'));
% files = files(~contains({files.name}, '191210'));
%files = files(~contains({files.name}, 'stim'));
%files = files(~contains({files.name}, '_4'));       % remove ones already done

batch_getRaster_fcn(files,refPeriod_ms);

%% plot before and after ttx
%match slice no. and spike detection method
fprintf(strcat('\n','\n','plotting before vs after TTX...','\n','\n'))

clearvars -except refPeriod_ms; close all

if ~exist('fs')
    fs = 25000;
end

bef_files = dir('*200127*mSpikes_4*.mat');                   % where your .mat files are
bef_files = bef_files(~~contains({bef_files.name}, 'Spikes')); %must contain 'Spikes')
bef_files = bef_files(~contains({bef_files.name}, 'TTX'));
%bef_files = bef_files(~contains({bef_files.name}, 'ttx'));
%bef_files = bef_files(~contains({bef_files.name}, 'Slice9'));
%bef_files = bef_files(~contains({bef_files.name}, 'Slice10'));
%bef_files = bef_files(~contains({bef_files.name}, 'Slice11'));
%bef_files = bef_files(~contains({bef_files.name}, 'Slice12'));
%bef_files = bef_files(~contains({bef_files.name}, 'Slice13'));
%bef_files = bef_files(~contains({bef_files.name}, 'Slice14'));
%bef_files = bef_files(~contains({bef_files.name}, 'Slice15'));
%bef_files = bef_files(~contains({bef_files.name}, 'Slice16'));
%bef_files = bef_files(~contains({bef_files.name}, 'Slice17'));
%bef_files = bef_files(~contains({bef_files.name}, 'adjM'));
%bef_files = bef_files(11); %for testing on one file; comment out for group analysis

for i = 1:length(bef_files)
    %beforeLabel     =    'before TTX';
    %afterLabel      =    'after TTX';
    %files = dir('*FTD*Group*Spikes*.mat*');
    %files = files(~contains({files.name}, 'slice1'));
    load(bef_files(i).name)
    
    if ~~exist('aSpikes')
        spikeMatrix = full(aSpikes);
        
    elseif ~~exist('cSpikes')
        spikeMatrix = full(cSpikes);
        
    elseif ~~exist('mSpikes')
        spikeMatrix = full(mSpikes);
        
    else
        disp('error - cant get spike mat')
    end
    clear aSpikes cSpikes mSpikes
    
    %stats before ttx
    stats.sp_count          = sum(sum(spikeMatrix));
    stats.num_active        = length(find(sum(spikeMatrix)>=1));
    
    scatter.mean_rate       = mean(sum(spikeMatrix))/(length(spikeMatrix)/fs);
    scatter.sem_rate        = std(sum(spikeMatrix))/sqrt(59);
    scatter.points.rate     = sum(spikeMatrix);
    
    before_length = length(spikeMatrix);
    
    title_strings = ...
        {{'Total number of spikes across all electrodes','(controlled for differences in recording length)'};...
        'Number of active electrodes'...
        };
    
    %get after file/s
    fileName = bef_files(i).name;
    %find file that contains ttx/TTX AND fileName
    spike_suffix_index = strfind(fileName,'Spikes')-3; %-3 to remove _aS in aSpikes for example
    spike_suffix = fileName(spike_suffix_index+1:end-4);
    aft_files = dir(strcat(fileName(1:spike_suffix_index),'*TTX*',spike_suffix,'.mat'));
    
    for j = 1 : length(aft_files)
        load(aft_files(j).name)
        if      ~~exist('aSpikes')
            spikeMatrix = full(aSpikes);
            
        elseif  ~~exist('cSpikes')
            spikeMatrix = full(cSpikes);
            
        elseif ~~exist('mSpikes')
            spikeMatrix = full(mSpikes);
            
        else
            disp('error - cant get spike mat')
        end
        clear aSpikes cSpikes mSpikes
        %stats after TTX
        %control for different length of recording
        after_length = length(spikeMatrix);
        scalar = before_length/after_length;
        
        stats.sp_count(j+1)             = sum(sum(spikeMatrix)) * scalar;
        stats.num_active(j+1)           = length(find(sum(spikeMatrix)>=1));
        
        %scatter.mean_rate(j+1)          = mean(sum(spikeMatrix))/(length(spikeMatrix)/fs);
        %scatter.sem_rate(j+1)           = std(sum(spikeMatrix))/sqrt(59);
        %scatter.points.rate(j+1)        = sum(spikeMatrix);
        %59 elecs; need to come up with a way to compensate for grounded
        %elecs automatically %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    num_stats = numel(fieldnames(stats)); %how many variables / after files
    %plot the figures for each stat if there are after files
    if length(aft_files) >= 1
        %stat_cell = struct2cell(stats);
        x = ['before';repmat(['after '],length(aft_files),1)];
        
        %%%%%%%%%%%%%%%%%% bar charts %%%%%%%%%%%%%%%%%%%%
        for j=1:num_stats
            stat_cell = struct2cell(stats);
            y = stat_cell{j};
            
            figure
            bar(y)
            aesthetics
            xticklabels(x)
            set(gca,'fontsize',20)
            ylabel(title_strings{j}, 'FontSize',14)
            xlabel('TTX condition')
            %ylim([0 200])
            if  strfind(fileName,'_') %remove underscores for title
                fileName1=fileName(1:end-4);
                fileName1(strfind(fileName,'_'))=' ';
                fileName1=strcat('{',fileName1,'}');
                title({fileName1,[' ']},'FontSize',16);
            else
                title({fileName(1:end-4),[' ']},'FontSize',16);
            end
            
            %%%%%%%%%%%%%%%%%%
            
            %save figure
            fprintf(strcat('\n','\n',bef_files(i).name(1:end-4),' saving bar chart...', '\n','\n'))
            try
                saveas(gcf,strcat(fileName(1:end-4),'_bar_',title_strings{j},'.png'));
            catch
                %if the title has multiple lines just save with first line
                %if you have cell within a cell 'title_strings{j}' will return
                %another cell rather than a string which causes error
                saveas(gcf,strcat(fileName(1:end-4),'_bar_',title_strings{j}{1},'.png'));
            end
            close all;
        end
        
        %%%%%%%%%%%%%%%%%%%% scatter plot %%%%%%%%%%%%%%%%%%%%
        %{
        figure
        %define x and y
        x_vals = linspace(1,1+length(aft_files),1+length(aft_files));
        x_vals1 = x_vals .* repmat(linspace(0.1,0.59,59),1+length(aft_files),1)';
        
        data = bsxfun(@times, rand(5,3), [50 150 100]);                     % Create Data
        dmean = mean(data);                                                 % Mean
        dci  = std(data)*tinv(0.975,size(data,1)-1);                        % Confidence Intervals
        xt = [1:3];                                                         % X-Ticks
        xtd = repmat(xt, size(data,1), 1);                                  % X-Ticks For Data
        sb = [xt'-ones(size(data,2),1)*0.1,  xt'+ones(size(data,2),1)*0.1]; % Short Bar X
        lb = [xt'-ones(size(data,2),1)*0.2,  xt'+ones(size(data,2),1)*0.2]; % Long Bar X
        figure(1)
        plot(xt, data, '+')
        hold on
        for k1 = 1:size(data,2)
        plot(lb(k1,:), [1 1]*dmean(k1), sb(k1,:), [1 1]*(dmean(k1)-dci(k1)), sb(k1,:), [1 1]*(dmean(k1)+dci(k1)), '-k')
        end
        hold off
        set(gca, 'XTick', xt, 'XTickLabel', {'Left','Middle','Right'})
        xlabel('Group')
        ylabel('Velocity (Furlongs/Fortnight)')
        
        aesthetics
        xticklabels(x)
        set(gca,'fontsize',20)
        ylabel(title_strings{j}, 'FontSize',14)
        xlabel('TTX condition')
        %ylim([0 200])
        if  strfind(fileName,'_') %remove underscores for title
            fileName1=fileName(1:end-4);
            fileName1(strfind(fileName,'_'))=' ';
            fileName1=strcat('{',fileName1,'}');
            title({fileName1,[' ']},'FontSize',16);
        else
            title({fileName(1:end-4),[' ']},'FontSize',16);
        end
        
        %%%%%%%%%%%%%%%%%%
        
        %save figure
        fprintf(strcat('\n','\n',bef_files(i).name(1:end-4),' saving scatter...', '\n','\n'))
        try
        saveas(gcf,strcat(fileName(1:end-4),'_scatter_','firing_rate_msem','.png'));
        catch
            %if the title has multiple lines just save with first line
            %if you have cell within a cell 'title_strings{j}' will return
            %another cell rather than a string which causes error
        saveas(gcf,strcat(fileName(1:end-4),'_bar_','firing_rate_msem','.png'));
        end
        %}
    end
    close all;
    clear fileName fileName1
    
end

%beforeFile      =
%afterFile       =
%bar_chart_organoid_basic_fcn(beforeFile,afterFile,beforeLabel,afterLabel)

%% get adjacency matrices
close all; clearvars -except refPeriod_ms

files = dir('*slice*mSpikes_3.mat');                   % where your .mat files are
files = files(~contains({files.name}, 'TTX'));
files = files(~contains({files.name}, 'ttx'));
files = files(~contains({files.name}, 'stim'));
files = files(~~contains({files.name}, 'FTD'));
files = files(~contains({files.name}, 'adjM'));
files = files(6:end);
%%%%%%% set these parameters first:
method = 'tileCoef';
sync_win_s = 0.200; %synchroncity window in seconds; e.g. 1 is +/- 1 s is considered synchronous (2DeltaT = 2 s)
rec_length_s = 360;
fs = 25000;
rec_length_samps = fs * rec_length_s;

%%%% downsampling:
num_samples = 1000 * rec_length_s; %defaults to 1 sample per ms
ds_rate = num_samples/rec_length_samps; %scalar for time_window %downsampling rate
sync_win = sync_win_s * ds_rate; %downsample the sync_win


fprintf(strcat('\n','\n','getting adjacency matrices...','\n','\n'))

% batch_getAdj_fcn(method,files,sync_win,num_samples,ds_rate);
% if recording lengths are all the same, uncomment above line and set
% rel_length_s above; otherwise, use below code to automatically find
% recording lengths.

for i = 1:length(files)
    files1 = files(i);
    spikes = struct2cell(load(files(i).name,'*Spikes'));
    rec_length_s = length(spikes{1}) / fs;
    rec_length_samps = fs * rec_length_s;
    num_samples = 1000 * rec_length_s; %defaults to 1 sample per ms
    ds_rate = num_samples/rec_length_samps; %scalar for time_window %downsampling rate
    sync_win = sync_win_s * ds_rate; %downsample the sync_win
    
    batch_getAdj_fcn(method,files1,sync_win,num_samples,ds_rate);
end

%% load adjMs and save as PNG
clear files adjM
files = dir('*slice*mSpikes_3_*adjM_0.175*.mat'); 
files = files(~~contains({files.name}, 'FTD'));
files = files(~contains({files.name}, 'CTRL'));
files = files([16 22 28 30]);

for i=1:length(files)
    fileName    =   files(i).name;
    if ~exist(strcat(fileName(1:end-4),'_connectivity_matrix_',num2str(sync_win/ds_rate),'.png'))
        load(fileName);
        if ~~exist('adjM2') & ~exist('adjM')
            adjM = adjM2;
        end
        figure; imagesc(adjM);
        aesthetics
        ylabel('Electrode')
        xlabel('Electrode')
        cb = colorbar;
        % ylabel(cb, 'Spike count')
        ylabel(cb, 'Correlation')
        cb.TickDirection = 'out';
        % cb.Ticks = 0:5; % for slice 5 specifically
        set(gca,'TickDir','out');
        %cb.Location = 'Southoutside';
        cb.Box = 'off';
        set(gca, 'FontSize', 14)
        caxis([0,1])
        %determine file name and title of fig
        fileName1=files(i).name;
        if  contains(files(i).name,'_') %remove underscores for title
            fileName1(strfind(fileName1,'_'))=' ';
            %fileName1=strcat('{',fileName1,'}');
            title({strcat(fileName1(1:end-4),' connectivity matrix'),' '});
        else
            title({strcat(files(i).name(1:end-4),' connectivity matrix'),' '});
        end
        %make title font smaller
        ax = gca;
        ax.TitleFontSizeMultiplier = 0.7;
        
        % add white lines
        hold on
        for j = 1:length(channels)
            plot([j - 0.5 j - 0.5], [0.5 60.5],'LineStyle','-','LineWidth',0.5,'Color',0.7*[1 1 1])
            plot([0.5 60.5]       , [j - 0.5 j - 0.5],'LineStyle','-','LineWidth',0.5,'Color',0.7*[1 1 1])
        end
        hold off
        
        %save as PNG
        fprintf(strcat('\n','\n',files(i).name(1:end-4),' saving adjM...', '\n','\n'))
        saveas(gcf,strcat(fileName(1:end-4),'_connectivity_matrix.png'));
        close all; clear adjM fileName fileName1
    else
    end
end


%% save heatMaps
clearvars -except refPeriod_ms

% files = dir('191210_slice1_DIV_g04_2018*_2min*-0.250*.mat'); 
% files = dir('191210_slice1_DIV_g04_2018*_2min*s_3_*.mat'); 
files = dir('191210_slice1_DIV_g04_2018*_2min*aSp*.mat'); 
files = files(~~contains({files.name}, 'TTX'));

% files = files(1);

option = 'logr'; %logc, count or rate (gets capped around 10,000 spikes)
fprintf(strcat('\n','\n','getting heat maps of:','\n',option,...
    ' spikes' ,'\n','\n'))
batch_getHeatMaps_fcn(files,option)

%% identify most spikey channel and period within that channel
close all; clearvars -except refPeriod_ms
files = dir('*slice*mSpikes_3.mat'); 
files = files(~contains({files.name}, 'FTD'));

fprintf(strcat('\n','\n','plotting spike overlays and marked filtered traces',...
    '\n','for spikiest channel' ,'\n','\n'))

option = 'diagonal';
progressbar('files')
for i = 1:length(files)
    % get spike mat
    fileName = files(i).name;
    spike_suffix_index = strfind(fileName,'Spikes')-3; %-3 to remove _aS in aSpikes for example
    spike_suffix = fileName(spike_suffix_index+1:end-4);
    
    % adaption for in case filename contans ref period 
    if ~isempty(strfind(spike_suffix,'RF'))
        RFindex = strfind(spike_suffix,'RF')-2;
        spike_suffix = spike_suffix(1:RFindex);
    end
    
    if      strcmp(spike_suffix(2:3),'mS')
        method = 'Manuel';
        parameter = str2double(spike_suffix(10:end));
    elseif  strcmp(spike_suffix(2:3),'cS')
        method = 'cwt';
        parameter = str2double(spike_suffix(11:end));
    elseif  strcmp(spike_suffix(2:3),'aS')
        method = 'abs';
        parameter = str2double(spike_suffix(end-7:end));
    else
        disp('error! Cannot determine spike detection mehtod for overlay')
    end
    % below functoin currently ignores artefacts for plotting;
    % create script that removes artefacts before spike detection or
    % removes spikes during spike detection if amplitude suggests it's
    % artefactual
    spike_overlay_fcn(fileName,method,parameter,refPeriod_ms,option);
    progressbar(i/length(files))
    %need to add correction to this function for where there are 0 spikes
end

%% plot graph features
close all; clearvars -except refPeriod_ms
% files = dir('*FTD*slice*_mSpikes_3_adjM_0.175.mat');
% files = files(~contains({files.name}, '1908'));
% files = files([17]);
scriptsDir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
filesDir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\Mecp2';
files = dir([filesDir '\' '*200205_4A*0.06*adjM*.mat']);

addpath(genpath(scriptsDir));

edge_thresh = 0.6;
node_thresh = edge_thresh;
progressbar('graph figs')
for file = 1:length(files)
    filename = files(file).name;
    plot_degree_weight_fcn(filesDir, filename);
    plot_rc_bc_fcn(filesDir, filename);
    plot_degree_edgenetwork_fcn(filesDir, filename,edge_thresh,node_thresh);
    progressbar(file/length(files));
end

%% reorder adjMs in order to use network plot in R

% %currently obsolete due to new plotting function in Matlab

% re_order_adjMs = 0; % for now don't reorder; do this manually later
% 
% clear all;close all
% if re_order_adjMs == 1;
%     files = dir('*200114*Slice6*mSpikes_5*adjM*.mat*');  % where your .mat files are
%     files = files(~contains({files.name}, 'reord'));
%     reorder_adjM_fcn(files);
% else
% end

%% to add: part that calls R script, plots network then save as PNG

