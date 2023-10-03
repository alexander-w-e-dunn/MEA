%% Calculating energy and ks statistics
% Written by Danyal Akarca and Alex Dunn,  University of Cambridge

%% About this script

%{

This script calculates energy and ks statistics 

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
datapath = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data\GNM_Data';
cd(datapath);

% add the brain connectivity toolbox to the path
brainconnectivity = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts\Other_scripts\2019_03_03_BCT';
addpath(brainconnectivity);

% % load spatial aspects
% load('MEAparameterspace');
% load('MEAcoordinates');
% load('MEAeuclidean');
% disp('Spatial aspects loaded');
%  
% % load organioid H9
% load('MEAweightedH9');
% load('MEAbinaryH9');
% load('MEAculturelistH9');
% load('MEAgenerativedataH9');
% load('MEAobservedH9');
% disp('Organoid H9 data loaded');

% load HD MEA data
% load('HDMEA1449gendata')
% load('HDMEA1449observed')
% load('HDMEA1449parameterspace');
% load('HD_MEA_binary');
% load('HD_MEA_binary_single');
% load('HD_MEA_Euclidean_single');
% load('HD_MEA_Euclidean');
% disp('HD MEA data loaded');

load('HDMEA_parameterspace_10k.mat');       % same for all timepoints

% we do not have the same n. cells at each time point for a given culture
% calculate euc. dist. mat. for each time point using
% HDMEA_analysis_revised_2021_02.m and data from Google Drive
load('HD_MEA_Euclidean_timepoint5.mat');     
load('HD_MEA_Gen_data_timepoint_5.mat');    
load('HD_MEA_Obs_data_timepoint_5.mat');

disp('HD MEA data loaded');

%% Form networks, calculate energy and KS statistics

% set which generative data and observed variables to analyse
% MEAgenerativedata = MEAgenerativedataH9;
% MEAobserved = MEAobservedH9;
% MEAgenerativedata = HDMEA1449gendata; 
% MEAobserved = HDMEA1449observed;
% MEAparameterspace = HDMEA1449parameterspace;
% MEAeuclidean = HD_MEA_Euclidean_single;
MEAgenerativedata   = HD_MEA_Gen_data; 
MEAobserved         = HD_MEA_Obs_data;
MEAparameterspace   = parameter_space;
clearvars -except MEAgenerativedata MEAobserved MEAparameterspace HD_MEA_Euclidean

%% run calculations
% initialise
MEAenergy = zeros(length(MEAgenerativedata),13,length(MEAparameterspace));
MEAks     = zeros(length(MEAgenerativedata),13,length(MEAparameterspace),4);

% loop over cultures
for culture = 1:length(MEAgenerativedata);
    culturedata = MEAgenerativedata{culture};
    % get the euc. dist. mat. for the nodes/cells in this culture
    MEAeuclidean = HD_MEA_Euclidean{culture};
    % get number of nodes/cells in the network for this culture
    n = length(MEAobserved{culture}{1,1});
    % loop over rules
    disp(newline + "       --> culture " + num2str(culture) + newline); beep
    tic
    for rule    = 1:13;
        gendata = squeeze(culturedata(rule,:,:));
        nB      = size(gendata,2);
        % loop over each simulation
        for iB = 1:nB
            b = zeros(n);
            b(gendata(:,iB)) = 1;
            b = b + b';
            y = cell(4,1);
            y{1} = sum(b,2);
            y{2} = clustering_coef_bu(b);
            y{3} = betweenness_bin(b)';
            y{4} = MEAeuclidean(triu(b,1) > 0);
            % calculate KS statistic
            K = [];
            for j = 1:4
                K(j) = fcn_ks(MEAobserved{culture}{j},y{j});
            end
            % save
            E = max(K,[],2);
            MEAenergy(culture,rule,iB) = E;
            MEAks(culture,rule,iB,:) = K;
            % display
            disp(sprintf('Culture %g Modeltype %g Simulation %g complete',culture,rule,iB));
        end
    end
    toc
end

%% save
% HDMEA1449energy     = MEAenergy;
% HDMEA1449ks         = MEAks;
% save('HDMEA1449energy.mat','HDMEA1449energy','-v7.3')
% save('HDMEA1449ks.mat','HDMEA1449ks','-v7.3')
timepoint = 5;
HD_MEA_energy     = MEAenergy;
HD_MEA_ks         = MEAks;
save(['HD_MEA_energy_timepoint_'    num2str(timepoint) '.mat'] , 'HD_MEA_energy','-v7.3')
save(['HD_MEA_ks_timepoint_'        num2str(timepoint) '.mat'] , 'HD_MEA_ks','-v7.3')

%% KS statistics function

function kstat = fcn_ks(x1,x2)
binEdges    =  [-inf ; sort([x1;x2]) ; inf];

binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
kstat = max(deltaCDF);
end