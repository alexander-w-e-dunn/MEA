%% Calulating distances between nodes: Generative models of across-scale neurodevelopment
% Written by Danyal Akarca and Alex Dunn,  University of Cambridge

%% About this script

%{
This script calculates the euclidean distances between nodes and visualises
them. The euclidean distance matrix will be an input to the generative
network models later.
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
datapath = '/Users/da04/Desktop/PhD/Research/GNM MEA/MEA data/';
cd(datapath);

% load MEA data
load('MEAweighted');
load('MEAbinary');
load('MEAculturelist');
load('MEAcoordinates');
load('MEAeuclidean');

%% 1. Visualise MEA distances

% node-wise summed distance
MEAeuclidean_sum                = sum(MEAeuclidean);
[sorted_MEAeuclidean_sum index] = sort(MEAeuclidean_sum);

% nine groups of nodes
group  = [1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9];
ugroup = unique(group);
% colours
colours   = flip(jet(length(ugroup)));

% plot
figure;
set(gcf, 'Position',  [1000, 1000, 1600, 400])
sgtitle('Multi-electrode array distances');
% plot euclidean distances
subplot(1,3,1); 
imagesc(MEAeuclidean);
xlabel('Sensor');
ylabel('Sensor');
c = colorbar;
c.Label.String = 'Euclidean distance'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';
set(gca,'TickLength',[0 0]);

% plot summed distance from each node
subplot(1,3,2); 
for k = 1:length(group);
    h = bar(k,sorted_MEAeuclidean_sum(k)); 
    h.FaceColor = colours(group(k),:); 
    hold on;
end
ylabel('Node-wise euclidean distance');
xticks([10:10:60]);
xlabel('Sensor');
set(gca,'TickLength',[0 0]);

% plot MEA
subplot(1,3,3);
for g = 1:length(group)
    scatter(MEAcoordinates(index(g),1),MEAcoordinates(index(g),2),80,colours(group(g),:),'filled'); 
    hold on;
end
xticks([1:8]); xlim([1 8]);
yticks([1:8]); ylim([1 8]);
xlabel('Sensor');
ylabel('Sensor');
set(gca,'TickLength',[0 0]);
