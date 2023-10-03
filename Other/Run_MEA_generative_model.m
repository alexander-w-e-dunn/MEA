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
datapath = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output\GNM_Data';
cd(datapath);

% add the brain connectivity toolbox to the path
brainconnectivity = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output\2019_03_03_BCT';
addpath(brainconnectivity);

% load spatial aspects
load('MEAparameterspace');
load('MEAcoordinates');
load('MEAeuclidean');

% load 14 day MEA WT
load('MEAweighted14');
load('MEAbinary14');
load('MEAculturelist14');

% load 14 day MEA HE
load('MEAweighted14HE');
load('MEAbinary14HE');
load('MEAculturelist14HE');

% load 21 day MEA WT
load('MEAweighted21');
load('MEAbinary21');
load('MEAculturelist21');

% load 28 day MEA WT
load('MEAweighted28');
load('MEAbinary28');
load('MEAculturelist28');

% load 28 day MEA KO
load('MEAweighted28KO');
load('MEAbinary28KO');
load('MEAculturelist28KO');

% load organoid EpiC (main WT)
load('MEAweightedEpiC');
load('MEAbinaryEpiC');
load('MEAculturelistEpiC');

% load organoid deficient WT
load('MEAweightedWTsli42');
load('MEAbinaryWTsli42');
load('MEAculturelistWTsli42');

% load organoid CS30 (FTD)
load('MEAweightedCS30');
load('MEAbinaryCS30');
load('MEAculturelistCS30');

% load organoid H9, removed culture 8 from the original
load('MEAweightedH9');
MEAweightedH9(:,:,8) = [];
load('MEAbinaryH9');
MEAbinaryH9(:,:,8) = [];
load('MEAculturelistH9');
MEAculturelistH9(8) = [];

%% 1. Set parameter space

% set vectors
vectoreta = linspace(-7,7,100);
vectorgam = linspace(-7,7,100);

% calculate grid
[a,b] = meshgrid(vectoreta,vectorgam);
c     = cat(2,a',b'); 

% calculate parameter space
parameter_space = reshape(c,[],2);

%% 2. Set generative rule hyperparameters

% *** set target networks ***
Atgt_set    = MEAbinary28;

% euclidean distance
D           = MEAeuclidean;

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
b = cell(size(Atgt_set,3),1);

% for KS calculations
observed = cell(size(Atgt_set,3),1);

for culture = 1:size(Atgt_set,3);
    
    % set target
    Atgt = squeeze(Atgt_set(:,:,culture));
    
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
  
    for i = 1:length(modeltype);
    
        % start clock
        tic;
        % run generative_model
        edges = generative_model(A,D,m,modeltype(i),modelvar,params);
        % keep value
        b{culture}(i,:,:) = edges;
        % stop clock
        t = toc;
        % display
        disp(sprintf('Culture %g, Modeltype %s complete (%s minutes)',culture,modeltype(i),t/60));
    
    end
    
    
end
