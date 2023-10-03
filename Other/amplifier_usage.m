% get amplifier usage
clear all

%find data directory in dropbox
cd 'C:\Users\alexd\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Mecp2 Project Jan 2019-'

% get mecp2 files
folders = dir('MPT*');

% for each folder containing a different prep get the recording files
for folder = 1 :length(folders)
    cd(folders(folder).name);
    files = dir('*.mcd');
    num_recs(folder) = length(files);
%     files = files(~[files.isdir]);
%     [~,idx] = sort([files.datenum]);
%     files = files(flip(idx));
%     [~,index] = sortrows({files.date}.'); files = files(index(end:-1:1)); clear index
    
    cd .. 
    clear files
end
num_recs = sum(num_recs);

% get the number of recording logs as I do one log for each session
% these are in their own folder
clearvars -except num_recs
cd 'C:\Users\alexd\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Mecp2 Project Jan 2019-\alex.recording.logs'
files = dir('*.txt');
% ignore files that aren't recording logs (I sometimes make notes on
% feeding etc)
files = files(~contains({files.name}, 'feeding','IgnoreCase',true));
files = files(~contains({files.name}, 'prep','IgnoreCase',true));

num_sessions = length(files)

% now repeat the above and get num reocrdings for other experiments
% (PV-Ai32 etc)
clearvars -except num_sessions num_recs
cd ..
folders = dir('PA*');

for folder = 1 :length(folders)
    cd(folders(folder).name);
    files = dir('*.mcd');
    num_recs2(folder) = length(files);
%     files = files(~[files.isdir]);
%     [~,idx] = sort([files.datenum]);
%     files = files(flip(idx));
%     [~,index] = sortrows({files.date}.'); files = files(index(end:-1:1)); clear index
    
    cd .. 
    clear files
end
num_recs2 = sum(num_recs2);

clearvars -except num_sessions num_recs num_recs2
cd ..
folders = dir('SMP*');

for folder = 1 :length(folders)
    cd(folders(folder).name);
    files = dir('*.mcd');
    num_recs3(folder) = length(files);
%     files = files(~[files.isdir]);
%     [~,idx] = sort([files.datenum]);
%     files = files(flip(idx));
%     [~,index] = sortrows({files.date}.'); files = files(index(end:-1:1)); clear index
    
    cd .. 
    clear files
end
num_recs3 = sum(num_recs2);

% sum the num recs for the different experiments
num_recs = num_recs + num_recs2 + num_recs3









