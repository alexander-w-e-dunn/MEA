%% Running generative models
% Written by Danyal Akarca and Alex Dunn,  University of Cambridge

%% About this script

%{

This script runs generative models for multi-electrode array data.

%}

%% 0. set paths and load required data

% clear workspace and command window
clear; clc;

% set latex formatting 
set(0,'DefaultTextInterpreter','Latex',...
    'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',1,...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

set(groot,'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');  

% set the current directory to the data
datapath = '/mhome/damtp/s/awed2/Documents/MatData/HD_MEA_Data';
cd(datapath);

% add the brain connectivity toolbox to the path
brainconnectivity = '/home/awed2/Documents/MatlabScripts/code_AD/2019_03_03_BCT';
addpath(brainconnectivity);

% % load spatial aspects
% load('MEAparameterspace');
% load('MEAcoordinates');
% load('MEAeuclidean');
% 
% % load 14 day MEA WT
% load('MEAweighted14');
% load('MEAbinary14');
% load('MEAculturelist14');
% 
% % load 14 day MEA HE
% load('MEAweighted14HE');
% load('MEAbinary14HE');
% load('MEAculturelist14HE');
% 
% % load 21 day MEA WT
% load('MEAweighted21');
% load('MEAbinary21');
% load('MEAculturelist21');
% 
% % load 28 day MEA WT
% load('MEAweighted28');
% load('MEAbinary28');
% load('MEAculturelist28');
% 
% % load 28 day MEA KO
% load('MEAweighted28KO');
% load('MEAbinary28KO');
% load('MEAculturelist28KO');
% 
% % load organoid EpiC (main WT)
% load('MEAweightedEpiC');
% load('MEAbinaryEpiC');
% load('MEAculturelistEpiC');
% 
% % load organoid deficient WT
% load('MEAweightedWTsli42');
% load('MEAbinaryWTsli42');
% load('MEAculturelistWTsli42');
% 
% % load organoid CS30 (FTD)
% load('MEAweightedCS30');
% load('MEAbinaryCS30');
% load('MEAculturelistCS30');
% 
% % load organoid H9, removed culture 8 from the original
% load('MEAweightedH9');
% MEAweightedH9(:,:,8) = [];
% load('MEAbinaryH9');
% MEAbinaryH9(:,:,8) = [];
% load('MEAculturelistH9');
% MEAculturelistH9(8) = [];

% load HD MEA data
% load('HD_MEA_binary_single');
load('HD_MEA_binary');
% load('HD_MEA_Euclidean_single');
load('HD_MEA_Euclidean');

%% 1. Set parameter space

% set vectors
vectoreta = linspace(-7,7,100);
vectorgam = linspace(-7,7,100);
% vectoreta = linspace(-7,7,10);
% vectorgam = linspace(-7,7,10);

% calculate grid
[a,b] = meshgrid(vectoreta,vectorgam);
c     = cat(2,a',b'); 

% calculate parameter space
parameter_space = reshape(c,[],2);

%% 2. Set generative rule hyperparameters

% *** set target networks ***
Atgt_set    = HD_MEA_binary; spike_sortd = 1;
% Atgt_set    = MEAbinary28; spike_sortd = 0;
%{
*** is data spike sorted ?
/ are there the same number of nodes in each network?
 if yes, put 1, otherwise put 0 ***

could add loop here: for Atgt_set = then list HD culture numbers
%}

% euclidean distance
% D           = MEAeuclidean;

% set the parameter for this subject
params      = parameter_space;

% calculate nparams
nparams     = length(params);
    
% model var
modelvar    = [{'powerlaw'},{'powerlaw'}];
    
% minimum edge
epsilon     = 1e-5;

% modeltype
modeltype = string({...
    'sptl',...          
    'neighbors',...     
    'matching',...      
    'clu-avg',...       
    'clu-min',...       
    'clu-max',...       
    'clu-diff',...      
    'clu-prod',...      
    'deg-avg',...       
    'deg-min',...       
    'deg-max',...       
    'deg-diff',...      
    'deg-prod'});      

%% Run generative models across the generative rules

% initialise

% for edge calculations
b = cell(length(Atgt_set),1);

% for KS calculations
observed = cell(length(Atgt_set),1);

for culture = 1:length(Atgt_set);
    
    % set target
    Atgt = Atgt_set{culture};
    
    % get distance matrix for spike_sorted data
    if spike_sortd == 1
        D = HD_MEA_Euclidean{culture}; % for multiple cultures
        % D = HD_MEA_Euclidean_single;
    end
    
    % calculate seed (idea for later, incorporate MST for some)
    ind = find(Atgt);
    k = randi(length(ind),1);
    randomedge = ind(k);
    A = zeros(size(Atgt));
    A(randomedge) = 1;
    A = A + A';
    
    % number of bi-directional connections
    m = nnz(Atgt)/2;
    
    % number of nodes
    n = length(Atgt);
    
    % calculate measures for later KS computations
    observed{culture}    = cell(4,1);
    observed{culture}{1} = sum(Atgt,2);
    observed{culture}{2} = clustering_coef_bu(Atgt);
    observed{culture}{3} = betweenness_bin(Atgt)';
    observed{culture}{4} = D(triu(Atgt,1) > 0);
    
    % initialise within culture
    b{culture} = zeros(length(modeltype),m,nparams);
  
    % if there are zeros not on the diagonal (i.e. cells in the same
    % location i.e. on top of each other, add small value
    if length(find(D == 0)) > length(D)
        % add 1e-05
        D = D + epsilon;
        % set diagonal to 0
        D( find(eye(size(D)) ==1) ) = 0;
    end
    
%     for i = 1:length(modeltype);
%     
%         % start clock
%         tic;
%         % run generative_model
%         edges = generative_model(A,D,m,modeltype(i),modelvar,params);
%         % keep value
%         b{culture}(i,:,:) = edges;
%         % stop clock
%         t = toc;
%         % display
%         disp(sprintf('Culture %g, Modeltype %s complete (%s minutes)',culture,modeltype(i),t/60));
%     
%     end
    
    
end

%% save

HD_MEA_generativedata_all = b;
save('HD_MEA_generativedata_all.mat','HD_MEA_generativedata_all','-v7.3')
save('HD_MEA_generativedata_all.mat','HD_MEA_generativedata_all','-v7.3')