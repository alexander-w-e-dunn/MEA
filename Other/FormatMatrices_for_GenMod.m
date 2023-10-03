%% load files and create 1 3d matrix
% Written by Alex Dunn and Danyal Akarca,  University of Cambridge

%{

This script: 

loads a list of .mat files containing MEA weighted adjacency matrices. Then
applies thresholding procedures and saves the different adjacency matrices,
weighted unthresholded, thresholded, thresholded and binarised all into one
3D matrix containing all MEA recordings. 

Note this currently works with up to 60x60 MEAs; this can increased by
setting the maximum matrix size. If the matrix is smaller than this, it
fills in the difference with 0s so that the matrices are all the same size.

%}
clear all
% load list of fileNames
% culture_list = readcell('GenModel_H9_101-619days_3DCulturesList.xlsx');
culture_list = readcell('GenModel_KO_2DCulturesList.xlsx');

matrix_suffix = '_DIV14_cSpikes_L-0.0627_RF1_adjM_0.05.mat';
% matrix_suffix = '_mSpikes_3_adjM_0.175.mat';

% SavefileName=strcat('MEA_matrices_H9_101-619days_human_iPSC');
SavefileName=strcat('MEA_matrices_KO_murine_DIV14');

for j=1  :  length(culture_list)
    load(strcat(culture_list{j,1},matrix_suffix));
    
    % added correction in case matrix is less than 60 x 60
    if length(adjM) ~= 60
        adjM(length(adjM) : 60 , length(adjM) : 60) = 0;
    end
    
    adjM_all(:,:,j)=adjM; 
    
    % 40% density proportion threshold
    adjM1 = adjM - 6 * eye(size(adjM));
    adjM1(find(isnan(adjM1))) = -5;
    w = find(adjM1 >-5);
    thr =   prctile(adjM1(w),60); % 40% density
    thr2 =  prctile(adjM1(w),75);
    thr3 =  prctile(adjM1(w),90);
    %         edge_thresh  =  mean([thr thr2 thr3]); % average across density (90th percentile is 10% density)
    edge_thresh = thr;  %threshold is 40% density; could floor to nearest 1 d.p. but loses accuracy (fig looks niceer though)
    % take first two decimal places
    edge_thresh = floor(edge_thresh * 100) / 100;
    if edge_thresh < 0.01
        edge_thresh = 0.01; % i.e. 1% above chance according to STTC; works better than 0
    end
    
    % could try absolute threshold of 0.3
%     edge_thresh = 0.3;
    
    binary_adjM = adjM1;
    binary_adjM(binary_adjM >= edge_thresh) = 1;
    binary_adjM(binary_adjM  < edge_thresh) = 0;
    
    binary_adjM_all(:,:,j) =  binary_adjM;
    
    include_nodes = find(sum(binary_adjM) >= 1);    
    trimmed_bin_adjM = binary_adjM(include_nodes,include_nodes); % remove nodes without edges
    trimmed_bin_adjM_all{j,1} = trimmed_bin_adjM;
    
    

    
    clear adjM binary_adjM adjM1 trimmed_bin_adjM
end

 
%% load channels and create 60 x 2 vector

 
    load(strcat(culture_list{1,1},matrix_suffix),'channels');
    coordstr=int2str(channels);
    for i=1:length(coordstr)
        
    coord(i,1)=str2double(coordstr(i));
    coord(i,2)=str2double(coordstr(i+60));
    end
    
%% save

save(SavefileName,'adjM_all','binary_adjM_all','binary_adjM_all','coord','culture_list','matrix_suffix','trimmed_bin_adjM_all')

% for all H9 ages together i manually extracted the ages from H9 org age
% list xls file and saved the vector below
% save(SavefileName,'ages','adjM_all','binary_adjM_all','binary_adjM_all','coord','culture_list','matrix_suffix','trimmed_bin_adjM_all')

