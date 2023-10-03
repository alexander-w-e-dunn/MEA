% analysis for MEA testing collaboration with Andras, Sagnik, Ilaria

%{

Note: exclude files containing stim artefacts and first run
remove_stim_artefacts_auto.m to get the cleaned files. Then you can run
this script on the cleaned files

%}

clear all
scriptsDir = 'C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
addpath(genpath(scriptsDir));
% dataDir = 'D:\MECP2_2019_AD\Scripts_and_Output\2.File_Conversion_Output';
rawDataDir = 'F:\Alex\MAT files\MEC';
spikeDataDir = 'C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\AnalysisOutput\By_culture\MEC\Spikes\cSpikes_L-0.0627_RP1';
refPeriod_ms = 1; %choose refractory period in ms

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
% clear all; close all
% DATA AND FILES
my_email = 'awed2@cam.ac.uk';
cd(rawDataDir)
% % TTX ONLY
% files = dir('MPT*.mat');
% files = files(~~contains({files.name}, 'TTX','IgnoreCase',true));
% files = files(~contains({files.name}, 'Spikes'));
% % files = files(~~contains({files.name}, 'AFTER'));
% files = files(3:4);

% DIV28 only:
files = dir('OWT*2202*.mat');
% files = files(~contains({files.name}, 'TTX','IgnoreCase',true));
files = files(~contains({files.name}, 'Spikes'));
% files = files(~contains({files.name}, 'DIV12'));
% files = files(~contains({files.name}, '200220'));
% % files = files(~contains({files.name}, 'DIV28'));
% files = files(~contains({files.name}, 'DIV4'));
% files = files(13);

% files = files(~contains({files.name}, 'DIV4'));
% % files1 = dir('MPT200107*DIV35.mat');
% files2 = dir('MPT200108*DIV35.mat');
% files3 = dir('MPT200109*DIV35.mat');
% files4 = dir('MPT200115*DIV*.mat');
%     files4 = files4(~contains({files4.name}, 'DIV49'));
%         files4 = files4(~contains({files4.name}, 'DIV42'));
%             files4 = files4(~contains({files4.name}, 'DIV35'));
% files5 = dir('MPT200205*DIV*.mat');
%     files5 = files5(~contains({files5.name}, 'DIV14'));
%         files5 = files5(~contains({files5.name}, 'DIV28'));
%             files5 = files5(~contains({files5.name}, 'Spikes'));
% files6 = dir('MPT200209*DIV*.mat');


% files = files(~contains({files.name}, 'stim'));
% files = files(~contains({files.name}, 'edited'));
% files = files(~contains({files.name}, '2min'));
% files = files(~contains({files.name}, 'MPT1'));
% files = files(~contains({files.name}, 'DIV49'));
% files = files(~contains({files.name}, '200209'));
% files = files(~contains({files.name}, '200220'));
% %  files = files(~contains({files.name}, 'DIV3'));
%  files = files(~contains({files.name}, 'DIV4'));
% files = files(~contains({files.name}, 'DIV5'));
% 
% files = files(39:end);

% files2 = dir('*MPT200115*DIV49*.mat');  
% files2 = files2(~contains({files2.name}, 'TTX','IgnoreCase',true));
% files2 = files2(~contains({files2.name}, 'Spikes'));
% files3 = dir('*MPT200205*DIV28*.mat');  
% files3 = files3(~contains({files3.name}, 'TTX','IgnoreCase',true));
% files3 = files3(~contains({files3.name}, 'Spikes'));
% files = [files;files2;files3];

%SPIKE DETECTION METHODS AND PARAMETERS
%option one: two methods, three parameters for each
%meths   ={'Manuel';'cwt'};
%params  =[4,5,6;                                % Threshold
%    -0.1254,0,0.1254]';                         % L parameter

%option two: two methods, one parameter for each
meths   ={'Manuel';'cwt'};
params  =[5,-0.0627];                               % threshold; L parameter

%option three: one method, one parameter
% meths   ={'Manuel'};                              % one method only for speed
% params  =[5,-0.0627];

%option four: one method, n parameters
% meths   ={'cwt'};
% params  =[4,4.25;                                % Threshold
%     0,-0.1880]';                       % L parameter
%
% CODE TO GET SPIKES AND SAVE:

fprintf(strcat('\n','\n','getting spike matrices...','\n','\n'))
disp(strcat({'number of files to do = '},num2str(length(files))))
progressbar('files','parameters','methods')

warning('off','MATLAB:load:variableNotFound') %turn warning off for now -
%just warns that a variable that i dont use anymore is not found

for m = 1:length(meths)
    method              =       meths{m}
    
    for p = 1:size(params,1);
        L               =       params(p+size(params,1))
        multiplier      =       params(p)
        
        batchGetSpike_function(rawDataDir,files,method,multiplier,L,refPeriod_ms,my_email)
        progressbar([],p/size(params,1),[]);       %update parameter progress
    end
    progressbar([],[],m/length(meths));            %update % methods done
end
warning('on','MATLAB:load:variableNotFound')

% move spike files to spikeDataDir
cd(rawDataDir)
movefile('*Spikes*',spikeDataDir)

% outputStats_organoid
% disp(' comment out line 25–27 and 89–90 from MEA test collab script!!')
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
%% identify most spikey channel and period within that channel
close all; clearvars -except refPeriod_ms spikeDataDir rawDataDir

cd(spikeDataDir)

files = dir('*cSpikes*RP1.mat');
% files = files(~~contains({files.name}, 'TTX','IgnoreCase',true));
% files = files(~~contains({files.name}, 'mSpikes'));
% files = files(~~contains({files.name}, 'BEFORE'));
files = files(~contains({files.name}, 'adjM'));
% files = files([12:13]);

% % DIV28 only
% files = dir('MPT*.mat');
% files = files(~contains({files.name}, 'TTX','IgnoreCase',true));
% files = files(~~contains({files.name}, 'Spikes'));
% files = files(~contains({files.name}, 'DIV12'));
% files = files(~~contains({files.name}, 'RF'));
% files = files(~~contains({files.name}, 'DIV28'));
% 
% files = files(93:94);

% files = dir('*MPT2001*_cSpikes_L0.mat*');  % where your .mat files are
% files = files(~contains({files.name}, 'TTX'));
% files = files(~contains({files.name}  , 'ttx'));
% %files = files(~contains({files.name}, 'Slice1.'));
% %files = files(~contains({files.name}, 'Slice1_'));
% %files = files(~contains({files.name}, 'Slice2'));
% %files = files(~contains({files.name}, 'Slice3'));
% %files = files(~contains({files.name}, 'Slice4'));
% files = files(~contains({files.name}, 'DIV07'));
% files = files(~contains({files.name}, '191210'));
% %files = files(~~contains({files.name}, 'Spikes')); %must contain 'Spikes')

fprintf(strcat('\n','\n','plotting spike overlays and marked filtered traces',...
    '\n','for spikiest channel' ,'\n','\n'))

% option = 'diagonal';
option = 'spikiest';
% option = [14];
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
        parameter = str2double(spike_suffix(10:strfind(spike_suffix,'RP')-2));
    elseif  strcmp(spike_suffix(2:3),'cS')
        method = 'cwt';
        parameter = str2double(spike_suffix(11:strfind(spike_suffix,'RP')-2));
    elseif  strcmp(spike_suffix(2:3),'aS')
        method = 'abs';
        if isempty(strfind(spike_suffix,'_based_on'))
            parameter = str2double(spike_suffix(10:strfind(spike_suffix,'RP')-2));
        elseif ~isempty(strfind(spike_suffix,'_based_on'))
            parameter = 0;
        end
    else
        disp('error! Cannot determine spike detection mehtod for overlay')
    end
    % below functoin currently ignores artefacts for plotting;
    % create script that removes artefacts before spike detection or
    % removes spikes during spike detection if amplitude suggests it's
    % artefactual
    spike_overlay_fcn(fileName,method,parameter,refPeriod_ms,option,rawDataDir,spikeDataDir);
    %need to add correction to this function for where there are 0 spikes
    progressbar(i/length(files))
end

%% creater rasters and save
close all; clearvars -except refPeriod_ms spikeDataDir rawDataDir files

cd(spikeDataDir)

fprintf(strcat('\n','\n','creating raster plots...','\n','\n'))
cd(spikeDataDir)
%inputs:

% TTX ONLY
% files = dir('MPT*.mat');
% files = files(~~contains({files.name}, 'TTX','IgnoreCase',true));
% files = files(~~contains({files.name}, 'aSpikes'));
% files = files(13:15);

% DIV28 ONLY
% files = dir('*MPT*2206*DIV*cSpikes*.mat');
% files = files(~~contains({files.name}, 'TTX','IgnoreCase',true));
% files = files(~~contains({files.name}, 'Spikes'));
% files = files(~contains({files.name}, 'adjM'));
% files = files(~contains({files.name}, 'DIV28'));

% PAT ONLY
% files = dir('PAT*17*cSpikes_L0.mat');


% files = dir('*FTD*mSpikes_3.mat');                   % where your .mat files are
%files = files(~contains({files.name}, 'stim'));
%files = files(~contains({files.name}, '_4'));       % remove ones already done
% files = files(41:end);

% files = dir('MPT200115_1D*cSpikes_L-0.0627*1.mat');
% % files = dir('MPT190403_6B*aSpikes_m*.mat');
% % files = files(~~contains({files.name}, 'TTX','IgnoreCase',true));
% % files = files(~~contains({files.name}, 'Spikes'));
% files2 = dir('MPT200108_4B*cSpikes_L-0.0627*1.mat');  
% files3 = dir('MPT190403_6C*cSpikes_L-0.0627*1.mat');  
% files = [files;files2;files3];
% files = files(~contains({files.name}, 'DIV07'));
% files = files(~contains({files.name}, 'DIV4'));
% files = files(~contains({files.name}, 'DIV5'));

batch_getRaster_fcn(files);
%% save heatMaps
close all; clearvars -except refPeriod_ms spikeDataDir rawDataDir files
cd(spikeDataDir)
% files = dir('MPT190403_6B*cSpikes_L-0.0627_RF1.mat');
% files = dir('*2112*cSpikes_L-0.0627_RP0.5.mat');
% files = dir('MPT190403_6B*aSpikes_m*.mat');
% files = files(~~contains({files.name}, 'TTX','IgnoreCase',true));
% files = files(~~contains({files.name}, 'Spikes'));
% files2 = dir('MPT200108_4B*cSpikes_L-0.0627*1.mat');  
% files3 = dir('MPT190403_6C*cSpikes_L-0.0627*1.mat');  
% files = [files;files2;files3];
% files = files(~contains({files.name}, 'DIV07'));
% files = files(~contains({files.name}, 'DIV4'));
% files = files(~contains({files.name}, 'DIV5'));
% 
% 
% files = files (23);

option = 'logr'; %logc, logr, count or rate (gets capped around 10,000 spikes)
fprintf(strcat('\n','\n','getting heat maps of:','\n',option,...
    ' spikes' ,'\n','\n'))
batch_getHeatMaps_fcn(files,option)
beep

%% plot before and after ttx
%match slice no. and spike detection method
fprintf(strcat('\n','\n','plotting before vs after TTX...','\n','\n'))

clearvars -except refPeriod_ms; close all

if ~exist('fs')
    fs = 25000;
end

bef_files = dir('*1217*cSpikes_L-0.0*.mat');                   % where your .mat files are
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
    % aft_files variable requires filenames to be identical + TTX on the
    % end
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
close all; clearvars -except refPeriod_ms spikeDataDir rawDataDir files
cd(spikeDataDir)
% files = dir('MPT200205_6*_DIV21*cSpikes_L-0.0627_RF1.mat');
% files = dir('*2112*cSpikes_L-0.0*.mat');
% files = dir('MPT*cSpikes_L-0.0627_RF1.mat');  % where your .mat files are
files = files(~contains({files.name}, 'adjM'));
% files = files(~contains({files.name}, 'TTX','IgnoreCase',true));
% files = files(~contains({files.name}, 'DIV12'));
% files = files(~contains({files.name}, 'DIV28'));
% files = files(~contains({files.name}, 'DIV07'));

%%%%%%% set these parameters first:
method = 'tileCoef';
% sync_win_s = [0.005 0.01 0.02 0.05]; %synchroncity window in seconds; e.g. 1 is +/- 1 s is considered synchronous (2DeltaT = 2 s)
sync_win_s = 0.01;
fs = 25000;

for f = 1 : length(files)
    filedata        = load(files(f).name);
    filename        = files(f).name;
    str_idx         = strfind(lower(filename),'spikes');
    spikes_string   = filename(str_idx-1:str_idx+5);
    rec_length_s    = length(filedata.(spikes_string))/fs;
    rec_length_samps = fs * rec_length_s;
    
    %%%% downsampling:
    num_samples = 1000 * rec_length_s; %defaults to 1 sample per ms
    ds_rate = num_samples/rec_length_samps; %scalar for time_window %downsampling rate
    
    fprintf(strcat('\n','\n','getting adjacency matrices...','\n','\n'))
    for i = 1:length(sync_win_s)
        
        sync_win = sync_win_s(i) * ds_rate; %downsample the sync_win
        
        batch_getAdj_fcn(method,files(f),sync_win,num_samples,ds_rate);
        
        clear sync_win 
    end
    sprintf('files done: %g',f)
    clear filedata filename str_idx spikes_string rec_length_s rec_length_samps
end

%% load adjMs and save as PNG
clear files adjM
files = dir('*M*adjM*.mat');
% files = dir('MPT200205_6*_DIV21*cSpikes_L-0.0627*adjM_0.05*.mat');  % where your .mat files are
% files = files(~contains({files.name}, 'CTRL'));
% files = files(~contains({files.name}, 'DIV28'));
% files = files(~contains({files.name}, 'TTX','IgnoreCase',true));
files = files(~~contains({files.name}, 'adjM'));
% files2 = dir('*MPT200115*DIV49*cSpikes_L0_adjM*.mat');  
% files3 = dir('*MPT200205*DIV28*cSpikes_L0_adjM*.mat');  
% files = [files;files2;files3];
% files = files(~contains({files.name}, 'reord'));

for i=1:length(files)
    fileName    =   files(i).name;
    dt_idx      = strfind(fileName,'_adjM_')+6;
    if ~exist(strcat(fileName(1:end-4),'_connectivity_matrix_',fileName(dt_idx:end-4),'.png'))
        load(fileName);
        try
            [~, I] = fcn_sort_channels_MEA(channels,'columns');
            adjM = adjM(I,I);
        catch
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
        if  contains(fileName1,'_') %remove underscores for title
            fileName1(strfind(fileName1,'_'))=' ';
            %fileName1=strcat('{',fileName1,'}');
            title({strcat(fileName1(1:end-4),' connectivity matrix'),' '});
        else
            title({strcat(fileName1(1:end-4),' connectivity matrix'),' '});
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
        saveas(gcf,strcat(fileName(1:end-4),'_connectivity_matrix.png')); %actual_sync_win=num2str(sync_win*(rec_length_samps/num_samples))
        close all; clear adjM fileName fileName1
    else
    end
end
beep

%% plot graph features
close all; clearvars -except refPeriod_ms files
% files = dir('MPT200220*_cSpikes_L-0.0627_RF1_adjM_0.05.mat');
% files = files(~~contains({files.name}, 'MPT190403_6B'));
% files = dir('MPT_220603_6C_DIV21*cSpikes_L-0.0*_adjM_*.mat');
% files = dir('MPT200205_6B*cSpikes_L-0.0627*1_adjM_0.05.mat');
% files = dir('MPT190403_6B*aSpikes_m*.mat');
% files = files(~~contains({files.name}, 'TTX','IgnoreCase',true));
% files = files(~~contains({files.name}, 'Spikes'));
% files2 = dir('MPT200108_4B*cSpikes_L-0.0627*1_adjM_0.05.mat'); 
% files2 = dir('MPT200115_2B*cSpikes_L-0.0627*1_adjM_0.05.mat');  
% % files3 = dir('MPT190403_6C*cSpikes_L-0.0627*1_adjM_0.05.mat');  
% files3 = dir('MPT190403_6B*cSpikes_L-0.0627*1_adjM_0.05.mat');  

% files = [files;files2;files3];
% files = files(~contains({files.name}, 'DIV07'));
% files = files(~contains({files.name}, 'DIV4'));
% files = files(~contains({files.name}, 'DIV5'));

directory = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\AnalysisOutput\By_culture\MPT\Spikes\cSpikes_L-0.0627_RP1';
edge_thresh = 0.01;
node_thresh = edge_thresh;
scale = 'logr';
cblim = [0.01 100];
thr_type = 'absolute';
thresholds = [0.1 0.3 0.5];


% progressbar('graph figs')
% for file = 1:length(files)
%     filename = files(file).name;
%         load(files(file).name)
%   
%     plot_degree_weight_fcn(directory, filename);
%     plot_rc_bc_fcn(directory, filename);
%     plot_degree_edgenetwork_fcn(directory, filename,edge_thresh,node_thresh);
% %     plot_FR_edgenetwork_fcn(directory, filename, edge_thresh, node_thresh, scale, cblim);
%     
%     progressbar(file/length(files));
% end


% network plots with firing rates, with dynamic threshold
progressbar('graph figs')
for file = 1:length(files)
    filename = files(file).name;
    load(files(file).name)
    adjM1 = adjM - 6 * eye(size(adjM));
    adjM1(find(isnan(adjM1))) = -5;
    w = find(adjM1 >-5);
    if strcmp(thr_type,'proportional')
        thr =   prctile(adjM1(w),60); % 40% density
        thr2 =  prctile(adjM1(w),75);
        thr3 =  prctile(adjM1(w),90);
        %         edge_thresh  =  mean([thr thr2 thr3]); % average across density (90th percentile is 10% density)
        edge_thresh = thr;  %threshold is 40% density; could floor to nearest 1 d.p. but loses accuracy (fig looks niceer though)
        % take first two decimal places
        edge_thresh = floor(edge_thresh * 100) / 100;
    end
%     absolute thresh:
    % edge_thresh = 0.8;
    
    if edge_thresh < 0 
        edge_thresh = 0.01; % i.e. 1% above chance according to STTC; works better than 0
    end
    node_thresh = edge_thresh;
    clearvars -except cblim directory edge_thresh files ...
        node_thresh refPeriod_ms scale filename file ...
        thr_type thresholds
    %         % sanity checks
    %         adjM2 = adjM1;
    %         adjM2(find(adjM2 < thr)) = 0;
    %         adjM3 = adjM2;
    %         adjM2(find(adjM2 >= thr)) = 1;
    %         f1=figure; imagesc(adjM1);
    %         f2=figure; imagesc(adjM2);
    %         f3 = figure; histogram(adjM3(find(adjM3 >= thr)));
    %         f4 = figure; histogram(adjM1(w));
    
%         plot_degree_weight_fcn(directory, filename);
%         plot_rc_bc_fcn(directory, filename,thr_type,thresholds);
%         plot_degree_edgenetwork_fcn(directory, filename,edge_thresh,node_thresh);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    could add in that only active channels are plotted
    this would require loading the spike file, finding the channels with
    <10 spikes, then using this vector of inactive channels to set all the
    weights associated with these channels to 0
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_FR_edgenetwork_fcn(directory, filename, edge_thresh, node_thresh, scale, cblim);
%     plot_thr_network_fcn(filename,'proportional',0.4,'weighted')
%     plot_binary_approach(filename,'proportional',0.4)

%     plot_binary_approach(filename,'absolute',0.2)
%     plot_thr_network_fcn(filename,'absolute',0.1,'weighted')
    progressbar(file/length(files));
end


%% reorder adjMs in order to use network plot in R

% %obsolete due to new Matlab plotting function

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

%% load excel file of filenames and genotypes and plot adjM over DIVs
clear all;close all

set(groot,'defaultAxesFontName','Arial')
set(groot,'defaultAxesFontSize',22)
set(groot,'defaultAxesLineWidth',2)
set(groot,'defaultFigurePosition',[700   435   520   350])

table = readtable('files_genotypes2.csv');
% files = sortrows(files,'rec','ascend');
adjM_suffix = '_cSpikes_L-0.0627_RF1_adjM_0.05.mat';
num_groups = flip( unique(table.gt) );
timepoints = {'DIV14','DIV21','DIV28','DIV35'} ;
cultures = unique(table.culture);
gt_list = table.gt;

for i = 1  : length(num_groups)
    for j = 1 : length(cultures)
        counts(j,i) = strcmp(num_groups{i},gt_list{j});
    end
end

n = sum(counts) 
% f1 = figure; f1.Position = [478    87   896   879];
thrOption = 'proportional';
thr = 0.7; type = 'weighted';
for i = 1 : length(cultures)

    for j = 1 : length(timepoints)
        
        fileName = strcat(cultures{i},'_',timepoints{j},adjM_suffix);
%         subplot(length(cultures),length(gt_list),i-1 + j)
        plot_thr_network_fcn(fileName,thrOption,thr,type)
        
    end
    
end












