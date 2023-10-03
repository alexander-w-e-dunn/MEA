% find files in directory one that are not in directory 2
dir1 = 'D:\MC_Rack Data\Alex\MCD files\Converted to raw\OWT';
dir2 = 'F:\Alex\MCD files\Converted to raw\OWT';
ext  = 'mcd';

cd(dir1)
file_list_1 = dir(['*.' ext])
file_list_1 = {file_list_1.name}'

cd(dir2)
file_list_2 = dir(['*.' ext])
file_list_2 = {file_list_2.name}'

unique_file_list = unique([file_list_1 ; file_list_2]);

% gives logical indices
% in1not2  = 1-ismember(unique_file_list, file_list_2);
% in2not1  = 1-ismember(unique_file_list, file_list_1);

in1not2 = unique_file_list(logical(1-ismember(unique_file_list, file_list_1)))
in2not1 = unique_file_list(logical(1-ismember(unique_file_list, file_list_2)))