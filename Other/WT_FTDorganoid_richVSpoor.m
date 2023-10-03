
% rich vs poor analysis and plot
clear all; close all
directory = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
files = dir('*slice*.mat');                   % where your .mat files are
files = files(~contains({files.name}, 'TTX'));
files = files(~contains({files.name}, 'ttx'));
files = files(~~contains({files.name}, 'mSpikes_3_adjM_0.175.mat'));
files = files(~contains({files.name}, 'stim'));
files = files(~~contains({files.name}, 'FTD'));
files = files(~contains({files.name}, 'spine'));
% files = files(~contains({files.name}, '190830'));
files = files(~contains({files.name}, '190814'));

for file = 1:length(files)    
    filename = files(file).name;
    [pc(file,:), bc(file,:), dk(file,:), ne(file,:), num_nodes(file,:), num_edges(file,:)] = richVSpoor_fcn(directory, filename);
end

groupdata = readcell('FTD_organoid_groups'); 

%% get groups

for file = 1:length(files)
    filename = files(file).name;
    for i = 1:length(groupdata)
        if strcmp(groupdata{i,1}, filename(1:end-25) )
            group_vec(file) = groupdata{i,2};            
        end
    end
end

% continued in R; manually extract and stack
% cell2table(groupdata);

%% plot 
addpath(genpath('gramm'));
X1 = [group_vec group_vec];
clear X
for i = 1   :   length(X1)
    grps = unique(X1);
    X(i,1) = find(X1(i) == grps);
end
Y = reshape(bc,[1 size(bc,1)*size(bc,2)])';
richness = [ones(1,length(bc)) zeros(1,length(bc))]';
load('example_data.mat');

% g=gramm('x',Model_Year,'y',MPG,'color',Cylinders,'subset',Cylinders~=3 & Cylinders~=5)
clear g
g = gramm('x',richness,'y',Y,'color',X);
g.geom_jitter()
g.stat_violin()
g.draw()






%% rich vs poor analysis and plot for synthetic matrices
clear all; close all
directory = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
files = dir('*slice*mSpikes_3*_adjM_0.175*.mat');                   % where your .mat files are
files = files(~contains({files.name}, 'TTX'));
files = files(~contains({files.name}, 'ttx'));
files = files(~contains({files.name}, 'stim'));
files = files(~~contains({files.name}, 'FTD'));
files = files(~~contains({files.name}, 'CTRL'));
% files = files(~contains({files.name}, '190830'));
files = files(~contains({files.name}, '190814'));

for file = 1:length(files)    
    filename = files(file).name;
    [pc(file,:), bc(file,:), dk(file,:), ne(file,:), num_nodes(file,:), num_edges(file,:)] = richVSpoor_fcn(directory, filename);
end

for file = 1:length(files)    
    filename = files(file).name;
%     plot_rc_bc_fcn(directory, filename);
    plot_degree_weight_fcn(directory, filename);
end

