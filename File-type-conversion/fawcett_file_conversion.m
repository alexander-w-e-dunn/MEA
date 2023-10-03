function fawcett_file_conversion(filename)
% directories
ext           = '.raw';
scriptsFolder = '/nfs/st01/hpc-cbi-sje30/awed2/scripts/Alex-PhD/file_conversion';
addpath(genpath(scriptsFolder));

inputFolder   = '/nfs/st01/hpc-cbi-sje30/awed2/file_conversion/raw_files';
outputFolder  = '/nfs/st01/hpc-cbi-sje30/awed2/file_conversion/mat_files';

% inputFolder   = 'F:\Alex\RAW files';
% outputFolder  = 'F:\Alex\MAT files';

cd(scriptsFolder);

%
MEA_load_bin_20210304(filename, 0, 'whole',inputFolder,outputFolder);
end