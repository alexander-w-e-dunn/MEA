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
% % load DIV14 WT
% load('MEAbinary14');
% load('MEAweighted14');
% load('MEAculturelist14');
% load('MEAgenerativedata14');
% load('MEAobserved14');
% disp('14 day wild type data loaded');
% 
% % load DIV14 HE
% load('MEAbinary14HE');
% load('MEAweighted14HE');
% load('MEAculturelist14HE');
% load('MEAgenerativedata14HE');
% load('MEAobserved14HE');
% disp('14 day heterozygous data loaded');
% 
% % load DIV21 WT
% load('MEAbinary21');
% load('MEAweighted21');
% load('MEAculturelist21');
% load('MEAgenerativedata21');
% load('MEAobserved21');
% disp('21 day wild type data loaded');
% 
% % load DIV28 WT
% load('MEAbinary28');
% load('MEAweighted28');
% load('MEAculturelist28');
% load('MEAgenerativedata28');
% load('MEAobserved28');
% disp('28 day wild type data loaded');
% 
% % load DIV28 KO
% load('MEAbinary28KO');
% load('MEAweighted28KO');
% load('MEAculturelist28KO');
% load('MEAgenerativedata28KO');
% load('MEAobserved28KO');
% disp('28 day knock-out data loaded');
% 
% % load organoid EpiC (main WT)
% load('MEAweightedEpiC');
% load('MEAbinaryEpiC');
% load('MEAculturelistEpiC');
% load('MEAgenerativedataEpiC');
% load('MEAobservedEpiC');
% disp('Organoid EpiC data loaded');
% 
% % load organoid EpiC (deficient WT)
% load('MEAweightedWTsli42');
% load('MEAbinaryWTsli42');
% load('MEAculturelistWTsli42');
% load('MEAgenerativedataWTsli42');
% load('MEAobservedWTsli42');
% disp('Organoid (deficient) EpiC data loaded');
% 
% % load organioid CS30 (FTD)
% load('MEAweightedCS30');
% load('MEAbinaryCS30');
% load('MEAculturelistCS30');
% load('MEAgenerativedataCS30');
% load('MEAobservedCS30');
% disp('Organoid CS30 data loaded');
% 
% % load organioid H9
% load('MEAweightedH9');
% load('MEAbinaryH9');
% load('MEAculturelistH9');
% load('MEAgenerativedataH9');
% load('MEAobservedH9');
% disp('Organoid H9 data loaded');

% load HD MEA data
load('HDMEA1449gendata')
load('HDMEA1449observed')
load('HDMEA1449parameterspace');
load('HD_MEA_binary');
load('HD_MEA_binary_single');
load('HD_MEA_Euclidean_single');
load('HD_MEA_Euclidean');
disp('HD MEA data loaded');

%% Form networks, calculate energy and KS statistics

% set which generative data and observed variables to analyse
% MEAgenerativedata = MEAgenerativedataH9;
% MEAobserved = MEAobservedH9;
MEAgenerativedata = HDMEA1449gendata; 
MEAobserved = HDMEA1449observed;
MEAparameterspace = HDMEA1449parameterspace;
MEAeuclidean = HD_MEA_Euclidean_single;

% set number of ROIs
% n         = 60;
n         = 252; % 252 cells in culture one of HD MEA data
%{ 
    need to save 'n' in generative data file for single HD MEA data and
    combined HD MEA data in order to have a variable number of ROIs/cells
%}

% initialise
MEAenergy = zeros(length(MEAgenerativedata),13,length(MEAparameterspace));
MEAks     = zeros(length(MEAgenerativedata),13,length(MEAparameterspace),4);

% loop over cultures
for culture = 1:length(MEAgenerativedata);
    culturedata = MEAgenerativedata{culture};
    % loop over rules
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
end

%% save
HDMEA1449energy     = MEAenergy;
HDMEA1449ks         = MEAks;
save('HDMEA1449energy.mat','HDMEA1449energy','-v7.3')
save('HDMEA1449ks.mat','HDMEA1449ks','-v7.3')

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