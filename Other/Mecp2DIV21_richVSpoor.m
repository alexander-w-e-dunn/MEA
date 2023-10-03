
% rich vs poor analysis and plot
clear all; close all
directory = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
list = readcell('files_included_gt.csv');
DIV = 35;
adjM_file_suffix = '_cSpikes_L-0.0627_RF1_adjM_0.05';
for i = 1 : length(list)-1
    files(i,1).name = strcat(list{i+1,1},'_DIV',num2str(DIV),adjM_file_suffix,'.mat');
end

for file = 1:length(files)    
    filename = files(file).name;
    [pc(file,:), bc(file,:), dk(file,:), ne(file,:), num_nodes(file,:), num_edges(file,:)] = richVSpoor_fcn(directory, filename);
end

groupdata = list(2:end,2); 

filenames = repmat(files,2,1);
num_nodes = repmat(num_nodes,2,1);
num_edges = repmat(num_edges,2,1);
bc(find(isnan(bc))) = 0;
pc(find(isnan(pc))) = 0;
ne(find(isnan(ne))) = 0;
dk(find(isnan(dk))) = 0;

for i = 1:length(filenames)
    
    RvsP{1,1} = 'filename';
    RvsP{1,2} = 'bc';
    RvsP{1,3} = 'pc';
    RvsP{1,4} = 'dk';
    RvsP{1,5} = 'ne';
    RvsP{1,6} = 'num_nodes';
    RvsP{1,7} = 'num_edges';
    
    RvsP{i+1,1} = filenames(i).name;
    RvsP{i+1,2} = bc(i);
    RvsP{i+1,3} = pc(i);
    RvsP{i+1,4} = dk(i);
    RvsP{i+1,5} = ne(i);
    RvsP{i+1,6} = num_nodes(i);
    RvsP{i+1,7} = num_edges(i);
    
end

RvsP(1,8) = {'richness'};
RvsP(2:length(files)+1       ,8) = {'rich'};
RvsP(2 + length(files) : end ,8) = {'poor'};

RvsP(1,9) = {'genotype'};
RvsP(2:length(files)+1       ,9) = groupdata(:);
RvsP(2 + length(files) : end ,9) = groupdata(:);


% Convert cell to a table and use first row as variable names
T = cell2table(RvsP(2:end,:),'VariableNames',RvsP(1,:));
 
% Write the table to a CSV file
writetable(T,strcat('RichVsPoor_MPT','_DIV',num2str(DIV),adjM_file_suffix,'.csv'))

% %% get groups
% 
% for file = 1:length(files)
%     filename = files(file).name;
%     for i = 1:length(groupdata)
%         if strcmp(groupdata{i,1}, filename(1:end-25) )
%             group_vec(file) = groupdata{i,2};            
%         end
%     end
% end
% 
% % continued in R; manually extract and stack
% % cell2table(groupdata);
% 
% %% plot 
% addpath(genpath('gramm'));
% X1 = [group_vec group_vec];
% clear X
% for i = 1   :   length(X1)
%     grps = unique(X1);
%     X(i,1) = find(X1(i) == grps);
% end
% Y = reshape(bc,[1 size(bc,1)*size(bc,2)])';
% richness = [ones(1,length(bc)) zeros(1,length(bc))]';
% load('example_data.mat');
% 
% % g=gramm('x',Model_Year,'y',MPG,'color',Cylinders,'subset',Cylinders~=3 & Cylinders~=5)
% clear g
% g = gramm('x',richness,'y',Y,'color',X);
% g.geom_jitter()
% g.stat_violin()
% g.draw()



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

