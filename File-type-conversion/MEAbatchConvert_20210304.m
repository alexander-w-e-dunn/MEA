function MEAbatchConvert_20210304(ext,inputFolder,outputFolder,scriptsFolder)
%form: MEA_batchConvert
%
%This function converts MEA bin (binary, RAW) files into mat files in a batch. It will
%attempt to convert all "raw" files in the directory, so make sure they are
%all MEA files or enter an extension to select a few.
%
%To create the RAW file, 
%1. open MC_DataTool
%2. File-Open Multiple
%3. Select files of interest
%4. Click "bin"
%5. Click "All"
%6. Make sure "Write header" and "Signed 16bit" are checked in lower right
%7. Click Save
%8. when done, click Close

% TODO: Improve the command line output of this

% Last update: 20180626 
% TS: Added the conversion options
% added directory inputs; AD 2021

%% AD - change directory into data folder...
% cd 'C:\Users\Windows\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Mecp2 Project Jan 2019-\RAW'
% cd 'C:\Users\Windows\Dropbox (Cambridge University)\NOG MEA Data\MEA Rich\WT\Left-Right'
cd(inputFolder)

%% Select conversion mode 


%convertOption = 'electrode'; % save electrode by electrode in a MEA-specific folder
 convertOption = 'whole'; % save the entire grid as one variable

%% initialize


if ~exist('ext','var')
    ext='.raw';
end;

ext2='.raw';

%% get files

d=dir;
files=[];

for i=1:length(d)
    if ~isempty(findstr(d(i).name,ext)) && ~isempty(findstr(d(i).name,ext2)) && length(d(i).name)>2 
        files=[files; i];
    end;
end;

files=d(files);

%% convert the files

plt=0; %AD added plt to be 1 so that plot is created; turn plot off during analysis for speed
numFails = 0;
for i=1:length(files)
    files(i).name
    skip=0;
    % find if file already converted (cd to output folder from
    % MEA_load_bin.m)
    % cd 'C:\Users\Windows\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Mecp2 Project Jan 2019-\MAT'
    % cd 'C:\Users\Windows\Dropbox (Cambridge University)\NOG MEA Data\MEA Rich\WT\Left-Right'
    cd(outputFolder)
    if ~~exist(strcat(files(i).name(1:end-4),'.mat'))
        %     for j=1:length(d)
        %         if ~isempty(strmatch(files(i).name(1:length(files(i).name)-4),d(j).name)) && d(j).isdir
        %             d(j).name
        skip=1;
        %     end;
    end;
    
    if skip==0
        % cd 'C:\Users\Windows\Dropbox (Cambridge University)\meaHeatMapOverlay\MatlabAnalysisScriptsAD'%cd to where scripts are
        % cd 'C:\Users\Windows\Dropbox (Cambridge University)\NOG MEA Data\MEA Rich\WT\Left-Right'
        cd(scriptsFolder)
        % add try loop: try convertion
        % if fails, catch and add the file name to a table of strings
        % print the list to a txt file and save as 'failed files'
        try
            MEA_load_bin_20210304(files(i).name, plt, convertOption,inputFolder,outputFolder); %AD added ' plt, ' as this is one of the required inputs
        catch
            numFails = numFails+1;
            fails{numFails,1} = files(i).name;
        end
    end
end;

if ~~exist('fails')
    formatOut = 'yyyy/mm/dd';
    datetoday = datestr(now,formatOut);
    datetoday = datetoday([1 2 3 4 6 7 9 10]);    
    writecell(fails,strcat(datetoday,'_files_failed_to_convert_',datestr(now, 'HH.MM'),'.txt'));
end

end

