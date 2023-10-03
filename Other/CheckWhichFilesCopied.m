clear all
cd 'C:\Users\alexd\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Mecp2 Project Jan 2019-\MAT files\Mecp2\Recordings'
files = dir('MPT*.mat');

cd 'D:\MECP2_2019_AD\Scripts_and_Output\AnalysisOutput\By_culture\MPT\Recordings'
files2 = dir('MPT*.mat');

files3 = [files; files2];

fnames1 = struct2cell(files);
fnames2 = struct2cell(files2);
% unique(fnames(1,:))

clear id ToCopy
ToCopy = {'Files to copy:'};
for file = 1 : length(fnames1)
    id{file} = find(contains(fnames2(1,:),fnames1{1,file}));
    if isempty(id{file})
       ToCopy{length(ToCopy)+1,1} =  fnames1{1,file};
    end
end