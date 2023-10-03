% takes data from cluster output from gen. model script and formats to be
% input to script to calculate energy. Creates a separate file for each
% timepoint. At each timepoint, the file will contain an n cultures x 1
% cell array
clear all
%% 
dataDir     = 'C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\GNM_Data\ForAlex\ForAlex\10pct\10pct\';
outputDir   = 'C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\GNM_Data';
timepoints  = 4:5;
cultures    = [1449; 1451; 1452; 1453; 1454; 1456; 1457];

for age = timepoints
    
    disp(newline + "       --> timepoint " + num2str(age) + newline)
    
    for culture = 1 : size(cultures,1)
        disp(newline + "       --> culture " + num2str(culture) + newline)

        cd([dataDir num2str(cultures(culture,:)) ])
        load(['edges-culture-' num2str(age)]);
        HD_MEA_Gen_data{culture,1} = table_edges;
        
        load( ['observed_' num2str(cultures(culture,:)) ] );
        HD_MEA_Obs_data{culture,1} = observed{age,1};
    end
    
    disp(newline + "       --> saving... " + newline)
    cd(outputDir);
    tic
    save(['HD_MEA_Gen_data_timepoint_' num2str(age)],'HD_MEA_Gen_data','-v7.3')
    save(['HD_MEA_Obs_data_timepoint_' num2str(age)],'HD_MEA_Obs_data','-v7.3')
    toc
    clear HD_MEA_Gen_data HD_MEA_Obs_data

end