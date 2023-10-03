%% 5. Longitudinal MEA energy landscapes: Generative models of across-scale neurodevelopment
% Written by Danyal Akarca and Alex Dunn,  University of Cambridge

%% About this script

%{

This script: 

1) Visualises longitudinal growth 

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

% add graphical functions that Alex has produced
graphpath = '/Users/da04/Desktop/PhD/Research/GNM MEA/MEA data/Alex functions/';
addpath(graphpath);

% set model types
modeltype = string({...
    'Spatial',...     % spatial model
    'Neighbors',...   % number of common neighbors
    'Matching',...    % matching index
    'C-Avg',...       % average clustering coeff
    'C-Min',...       % minimum clustering coeff
    'C-Max',...       % maximum clustering coeff
    'C-Diff',...      % difference in clustering coeff
    'C-Prod',...      % product of clustering coeff
    'D-Avg',...       % average degree
    'D-Min',...       % minimum degree
    'D-Max',...       % maximum degree
    'D-Diff',...      % difference in degree
    'D-Prod'})';      % product of degree

% set energy equation components
eqncomponents = string({...
    'Efficiency',...
    'Clustering',...
    'Degree',...
    'Edge length'});

% load spatial aspects
load('MEAparameterspace');
load('MEAcoordinates');
load('MEAeuclidean');

% load DIV14 WT
load('MEAbinary14');
load('MEAweighted14');
load('MEAculturelist14');
load('MEAgenerativedata14');
load('MEAobserved14');
load('MEAenergy14');
load('MEAks14');

% load DIV14 HE
load('MEAbinary14HE');
load('MEAweighted14HE');
load('MEAculturelist14HE');
load('MEAgenerativedata14HE');
load('MEAobserved14HE');
load('MEAenergy14HE');
load('MEAks14HE');

% load DIV21 WT
load('MEAbinary21');
load('MEAweighted21');
load('MEAculturelist21');
load('MEAgenerativedata21');
load('MEAobserved21');
load('MEAenergy21');
load('MEAks21');

% load DIV28 WT
load('MEAbinary28');
load('MEAweighted28');
load('MEAculturelist28');
load('MEAgenerativedata28');
load('MEAobserved28');
load('MEAenergy28');
load('MEAks28');

% load DIV28 KO
load('MEAbinary28KO');
load('MEAweighted28KO');
load('MEAculturelist28KO');
load('MEAgenerativedata28KO');
load('MEAobserved28KO');
load('MEAenergy28KO');
load('MEAks28KO');

%% 1. Set culture, model, form landscapes
 
culture = 7;
model   = 3;

% initialise
elandscape  = [];
kslandscape = []; 
energy      = {MEAenergy14 MEAenergy21 MEAenergy28};
ks          = {MEAks14 MEAks21 MEAks28};

% loop for each and keep
for div = 1:3
    % index into this div
    MEAenergy = energy{div};
    MEAks     = ks{div};
    % calculate landscapes
    se                 = squeeze(MEAenergy(culture,model,:));
    sp                 = MEAparameterspace;
    [u x y]            = unique(sp,'rows');
    se                 = se(x);
    se                 = reshape(se,[100 100]);
    % keep energy landscape
    elandscape(div,:,:) = se;
    for k = 1:4
        sk = squeeze(MEAks(culture,model,:,k));
        sk = sk(x);
        sk = reshape(sk,[100 100]);
        % keep ks landscape
        kslandscape(k,div,:,:) = sk;
    end
end

%% 2. Visualise network over time

% make graphs
g14 = graph(squeeze(MEAbinary14(:,:,culture)));
g21 = graph(squeeze(MEAbinary21(:,:,culture)));
g28 = graph(squeeze(MEAbinary28(:,:,culture)));

% number of edges
n14 = nnz(squeeze(MEAbinary14(:,:,culture)))/2;
n21 = nnz(squeeze(MEAbinary21(:,:,culture)))/2;
n28 = nnz(squeeze(MEAbinary28(:,:,culture)))/2;

% make landscapes
e14 = squeeze(elandscape(1,:,:));
e21 = squeeze(elandscape(2,:,:));
e28 = squeeze(elandscape(3,:,:));

% visualise
figure; 
set(gcf,'Position',[1000 1000 1200 600])
sgtitle(sprintf('MEA culture %g, %s',culture,modeltype(model)));

% cultures

% div14

subplot(2,3,1);

h14 = plot(g14,'EdgeAlpha',0.8,'EdgeColor','k','NodeColor','k','Xdata',MEAcoordinates(:,1),'Ydata',MEAcoordinates(:,2),'interpreter','latex');
title(sprintf('14 days, %g edges',n14));
xticks([]); yticks([]);
labelnode(h14,1:60,'');

subplot(2,3,4);
imagesc(e14);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
caxis([0 1]); xticks([]); yticks([]);
c = colorbar;
c.Label.String = 'Energy'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';

% div21

subplot(2,3,2); 

h21 = plot(g21,'EdgeAlpha',0.8,'EdgeColor','k','NodeColor','k','Xdata',MEAcoordinates(:,1),'Ydata',MEAcoordinates(:,2),'interpreter','latex');
title(sprintf('21 days, %g edges',n21));
xticks([]); yticks([]);
labelnode(h21,1:60,'');

subplot(2,3,5);
imagesc(e21);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
caxis([0 1]); xticks([]); yticks([]);
c = colorbar;
c.Label.String = 'Energy'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';

% div28

subplot(2,3,3); 

h28 = plot(g28,'EdgeAlpha',0.8,'EdgeColor','k','NodeColor','k','Xdata',MEAcoordinates(:,1),'Ydata',MEAcoordinates(:,2),'interpreter','latex');
title(sprintf('28 days, %g edges',n28));
xticks([]); yticks([]);
labelnode(h28,1:60,'');

subplot(2,3,6);
imagesc(e28);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
caxis([0 1]); xticks([]); yticks([]);
c = colorbar;
c.Label.String = 'Energy'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';

%% 3. Energy landscape changes

figure;
sgtitle('$$ \Delta $$E: Decreasing E corresponds to improved performance');
set(gcf,'Position',[1000 1000 1200 300])

% 14 to 21
subplot(1,3,1);
imagesc(e21-e14);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
xticks([]); yticks([]);
title('14 days to 21 days');
c = colorbar;
caxis([-0.5 0.5]);
c.Label.String = '$$ \Delta $$E'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';

subplot(1,3,2);
imagesc(e28-e21);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
xticks([]); yticks([]);
title('21 days to 28 days');
c = colorbar;
caxis([-0.5 0.5]);
c.Label.String = '$$ \Delta $$E'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';

subplot(1,3,3);
imagesc(e28-e14);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
xticks([]); yticks([]);
title('14 days to 28 days');
c = colorbar;
caxis([-0.5 0.5]);
c.Label.String = '$$ \Delta $$E'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';
