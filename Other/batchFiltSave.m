%script to run through and get filtered data from all files
% Author: Alex Dunn, Cambridge, 2019; last update: Nov. 2021
% Uses filtering approach in detectspikes.m  by Tim Sit, sitpakhang@gmail.com

clear all

%% load data
% %input data needs to be n samples x n channels
data_dir = 'G:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\By_culture\MPT\Recordings';
cd(data_dir);
save_dir = 'G:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\By_culture\MPT\Recordings\Developmental Dataset 2020';
% % get files containing this string
% files=dir('*190705*.mat');
% % remove files containing this phrase
% files = files(~contains({files.name}, 'Spikes'));
% files = files(~contains({files.name}, 'Filt'));
filelistcsv = 'files_genotypes_non_intense.csv';
data = readtable(filelistcsv);
% exclude data with intense periods
copyData = data(data.Intense_periods == 0,:);
files = copyData.Recording;

progressbar

%% Filter signal
for file=1:length(files)
    
    cd(data_dir);
    % check if input is struct array from 'dir' command
    % or cell array from csv file
    if iscell(files)
        raw_filename = [files{file},'.mat'];
    elseif isstruct(files)
        raw_filename = files(file).name;
    end
    
    % check if already filtered and saved, if not, load raw data
    if ~exist([raw_filename(1:end-4), '_Filtered', '.mat'])

        tic;
        load(raw_filename,'dat','channels','fs');
        disp(strcat('loaded: ',raw_filename))
        
        %filter
        lowpass = 600;
        highpass = 8000;
        wn = [lowpass highpass] / (fs / 2);
        filterOrder = 3; %changed from 3 to 5 by AD; and back
        [b, a] = butter(filterOrder, wn);
        filteredMatrix = filtfilt(b, a, double(dat));
        
        %save as .mat matrix n sample x n channel
        fileName = [raw_filename(1:end-4), '_Filtered', '.mat'];
        % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
        disp('saving...')
        cd(save_dir)
        save(fileName, 'filteredMatrix','channels','fs', '-v7.3');
        
        toc
        
        clear dat
        clear filteredMatrix
        clear channels
    else
        disp(strcat(files(file).name(1:end-4),'_file done already'))
    end
    progressbar(file/length(files));
end

