clear all
% transnetyx cost calculations
% dataDir = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Lab_and_supplies\Mouse management\Genotyping\Transnetyx and cage costs 2022-03-01 to 2022-08-31';
inputdata = inputdlg({['Paste directory here; e.g. C:\Users\alexd\Downloads' ,newline,...
    'This script will find use csv files with names starting with "OrderResults"'],
    ['Input spreadsheet name here e.g. Cage and Genotyping Costs Last 6 months.xlsx',...
    newline,'Please include file extension e.g. .xlsx']}); 
    
dataDir = inputdata{1};
spreadsheetname = inputdata{2};

% add \ to end of dataDir if needed
if ~strcmp(dataDir(end),'\') & ~strcmp(dataDir(end),'/') & contains(dataDir,'\')
    dataDir(end+1) = '\';
elseif ~strcmp(dataDir(end),'\') & ~strcmp(dataDir(end),'/') & contains(dataDir,'/')
    dataDir(end+1) = '/';
end

files = dir('OrderResults*.csv');
% input minimum cost per sample tested 
BasicCost   = 6.45;
% input additional cost per strain tested
AddCost     = 1.5;


cd(dataDir)

ORDER_NUMBER        = [];
STRAIN_NAME         = [];
NUM_SAMPLES         = [];
NUM_PROBES          = [];
PROBES_PER_SAMPLE   = [];
SAMPLES_OMITTED     = [];
COST                = [];

% set up an output table
% output 1 columns: order number; strain; probes tested; samplest tested; samples omitted; probes per sample; cost
% ouptut 2 columns: strain; number of probes; total samples tested across orders; total probes tested; total cost
% for output 2; get unique strains of output 1; loop through and get indices of that strain and sum the number of samples, probes and cost

for i = 1 : length(files)
    files(i).name
    % load data
    try
        data = readtable(files(i).name);
        % get well plates included in order (usually one but can have multiple)
        WellPlates = table2array(unique(data(:,'WellPlate')));
    catch
        data = readcell(files(i).name);
        data = cell2table(data(2:end,:),...
            "VariableNames",data(1,:))
        WellPlates = table2array(unique(data(:,'WellPlate')));
    end
    % get index of second hyphen in file name
    idx = strfind(files(i).name,'-'); idx=idx(2)-1;
    % get order number (14th character to character before first hyphen)
    OrderNumber = files(i).name(14:idx);
    % get list of strains
    strains = unique(data.Strain);
    % for each strain, get number of probes tested
    for j = 1 : length(strains)
        strain_idx = find(strcmp(data.Strain,strains{j}));
        % convert table to cell array
        C=table2cell(data(strain_idx,:));
        % remove first 5 columns that don't relate to different probes
        C = C(:,6:end);
        % convert all cells to strings
%         out=cellfun(@num2str,{C{:,6:end}},'un',0)
        % get number of probes tested as the number of results indicated by
        % + or -
%         NumProbes       = length([ strfind(strcat(out{:}),'+') strfind(strcat(out{:}),'-')] );
        NumProbes       = sum(~cellfun(@isempty,C) , 'all');
        % get number of tissue samples tested
        % number of samples in plate:
        NumSamples      = size(C,1);
        ProbesPerSample = ceil(NumProbes/NumSamples);
        % get number of samples and probes omitted and subtract this from costs
        % find which samples have all probes omitted (i.e. sample not
        % tested at all so doesn't pay the base cost
        SamplesOmitted  = sum(sum(strcmp(C,'Omit'),2) == ProbesPerSample);
        SamplesComplete = NumSamples - SamplesOmitted;
        ProbesOmitted   = sum(~~strcmp(C,'Omit'),'all');
        ProbesComplete  = NumProbes  - ProbesOmitted;
        
        totalCost   = BasicCost * SamplesComplete + (AddCost * ProbesComplete);
        
        ORDER_NUMBER        = [ORDER_NUMBER         ; string(OrderNumber)];
        STRAIN_NAME         = [STRAIN_NAME          ; string(strains{j})];
        NUM_SAMPLES         = [NUM_SAMPLES          ; NumSamples];
        NUM_PROBES          = [NUM_PROBES           ; NumProbes];
        PROBES_PER_SAMPLE   = [PROBES_PER_SAMPLE    ; ProbesPerSample];
        SAMPLES_OMITTED     = [SAMPLES_OMITTED      ; SamplesOmitted];
        COST                = [COST                 ; totalCost];
        
        % add results to table; 
    end
end

%% save to existing table

T = table(ORDER_NUMBER,STRAIN_NAME,NUM_SAMPLES,NUM_PROBES,PROBES_PER_SAMPLE,SAMPLES_OMITTED,COST)

% writetable(T,'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Lab_and_supplies\Mouse management\Genotyping\Transnetyx and cage costs 2022-03-01 to 2022-08-31\Cage and Genotyping Costs Last 6 months.xlsx',...
%     'Sheet','Transnetyx order costs')
writetable(T,[dataDir,spreadsheetname],'Sheet','Transnetyx order costs')

