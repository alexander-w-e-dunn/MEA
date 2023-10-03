%% Visualising MEA energy landscapes: Generative models of across-scale neurodevelopment
% Written by Danyal Akarca and Alex Dunn,  University of Cambridge

%% About this script

%{

This script 

1) Calculates the energy distributions for each of the generative rules.
2) Provides the energy landscape across this space for eta and gamma.

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

% display length of time
disp('Loading takes ~5 minutes');

% load spatial aspects
load('MEAparameterspace');
load('MEAcoordinates');
load('MEAeuclidean');
disp('Spatial aspects loaded');

% load DIV14 WT
load('MEAbinary14');
load('MEAweighted14');
load('MEAculturelist14');
load('MEAgenerativedata14');
load('MEAobserved14');
load('MEAenergy14');
load('MEAks14');
disp('14 day wild type data loaded');

% load DIV14 HE
load('MEAbinary14HE');
load('MEAweighted14HE');
load('MEAculturelist14HE');
load('MEAgenerativedata14HE');
load('MEAobserved14HE');
load('MEAenergy14HE');
load('MEAks14HE');
disp('14 day heterozygeous data loaded');

% load DIV21 WT
load('MEAbinary21');
load('MEAweighted21');
load('MEAculturelist21');
load('MEAgenerativedata21');
load('MEAobserved21');
load('MEAenergy21');
load('MEAks21');
disp('21 day wild type data loaded');

% load DIV28 WT
load('MEAbinary28');
load('MEAweighted28');
load('MEAculturelist28');
load('MEAgenerativedata28');
load('MEAobserved28');
load('MEAenergy28');
load('MEAks28');
disp('28 day wild type data loaded');

% load DIV28 KO
load('MEAbinary28KO');
load('MEAweighted28KO');
load('MEAculturelist28KO');
load('MEAgenerativedata28KO');
load('MEAobserved28KO');
load('MEAenergy28KO');
load('MEAks28KO');
disp('28 day knock-out data loaded');

% load organoid EpiC (main WT)
load('MEAweightedEpiC');
load('MEAbinaryEpiC');
load('MEAculturelistEpiC');
load('MEAgenerativedataEpiC');
load('MEAobservedEpiC');
load('MEAenergyEpiC');
load('MEAksEpiC');
load('MEAagesEpiC');
disp('Organoid EpiC data loaded');

% load organoid WTsli42
load('MEAweightedWTsli42');
load('MEAbinaryWTsli42');
load('MEAculturelistWTsli42');
load('MEAgenerativedataWTsli42');
load('MEAobservedWTsli42');
load('MEAenergyWTsli42');
load('MEAksWTsli42');
load('MEAagesWTsli42');
disp('Organoid WTsli42 data loaded');

% load organoid CS30 (FTD)
load('MEAweightedCS30');
load('MEAbinaryCS30');
load('MEAculturelistCS30');
load('MEAgenerativedataCS30');
load('MEAobservedCS30');
load('MEAenergyCS30');
load('MEAksCS30');
load('MEAagesCS30');
disp('Organoid CS30 (FTD) data loaded');

% load organoid H9
load('MEAweightedH9');
load('MEAbinaryH9');
load('MEAculturelistH9');
load('MEAgenerativedataH9');
load('MEAobservedH9');
load('MEAenergyH9');
load('MEAksH9');
load('MEAagesH9');
disp('Organoid H9 data loaded');

%% 1. Set which cultures to use

MEAenergy = MEAenergyH9;
MEAks     = MEAksH9;
MEAbinary = MEAbinaryH9;

%% 2. Remove specific cultures manually that are not connected or too small
 
% 14 WT - culture 15
% 14 HE - culture 9, culture 15
% 21 WT - none
% 28 WT - none
% 28 KO - culture 13

% *** set cultures to remove ***
remove    = [];

% recalculate
MEAenergy(remove,:,:) = []; 
MEAks(remove,:,:,:) = []; 
MEAbinary(:,:,remove) = []; 

% calculate number of edges in each culture and rank
nnzculture = [];
for n = 1:size(MEAbinary,3);
    A = squeeze(MEAbinary(:,:,n));
    nnzculture(n) = nnz(A)/2;
end
[~,nnzrank] = sort(nnzculture);

%% 3. Calculate minimum e, minimum ks, parameters

% energy
[e ei]    = min(MEAenergy,[],3);

% ks statistics
[ks ksi]  = (min(MEAks,[],3));
ks        = squeeze(ks);
ksi       = squeeze(ksi);

% parameters
p1        = MEAparameterspace(ei,1); p1 = reshape(p1,[size(MEAenergy,1) 13]);
p2        = MEAparameterspace(ei,2); p2 = reshape(p2,[size(MEAenergy,1) 13]);

%% 4. Visualise findings

% culture energy across rules
figure; 
set(gcf,'Position',[1000 1000 1200 600])
imagesc(e);
title('Cortical MEA cultures');
ylabel('Culture');
xlabel('Generative Rule');
xticks([1:13]);
xticklabels(modeltype);
c = colorbar;
caxis([0 1]);
c.Label.String = 'Energy'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';
c.Ticks = [0 1];
set(gca,'TickLength',[0 0]);

% culture KS statistics across rules
figure; 
set(gcf,'Position',[1000 1000 1200 600])
for k = 1:4
    subplot(2,2,k);
    imagesc(ks(:,:,k));
    title(eqncomponents(k));
    ylabel('Culture');
    xlabel('Generative Rule');
    set(gca,'TickLength',[0 0]);
    caxis([0 1]);
end

% mean best energy
rearrange = [3 2 4 6 7 8 5 9 11 12 13 10 1];

% Boxplot across the top network of the whole sample
figure;  
set(gcf,'Position',[1000 1000 1200 600])
h = boxplot(e(:,rearrange),'plotstyle','compact','color',[0.3 0.3 0.3]); 
title('Energy distribution of the best network for each generative rule'); 
ylabel('Energy'); 
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]);
xticklabels(modeltype(rearrange));

a = get(get(gca,'children'),'children'); 
bi = [53:65];
li = [66:79];

% set colours
c = ['r' 'r' 'g' 'g' 'g' 'g' 'g' 'b' 'b' 'b' 'b' 'b' 'k']; 
c = flip(c);
for i = 1:length(c)
    set(a(bi(i)),'Color',c(i));
    set(a(li(i)),'Color',c(i));
end
set(gca,'TickLength',[0 0])

%% 5. Tabulate findings across the sample

networktable = table(modeltype,mean(e)',std(e)',mean(p1)',std(p1)',mean(p2)',std(p2)','VariableNames',...
    {'Rule','Mean E','Std E','Mean P1','Std P1','Mean P2','Std P2'});
disp(networktable);

%% 6. Plot energy landscapes: across cultures

% *** pick model *** 
model = 5;

% initialise
elandscape  = [];
kslandscape = [];

% loop for each and keep
for culture = 1:size(MEAenergy,1);
    
    se                 = squeeze(MEAenergy(culture,model,:));
    sp                 = MEAparameterspace;
    [u x y]            = unique(sp,'rows');
    se                 = se(x);
    se                 = reshape(se,[100 100]);
    
    % keep energy landscape
    elandscape(culture,:,:) = se;
    
    for k = 1:4
        sk = squeeze(MEAks(culture,model,:,k));
        sk = sk(x);
        sk = reshape(sk,[100 100]);
        
        % keep ks landscape
        kslandscape(k,culture,:,:) = sk;
    end
    
end

% mean across participants
meanelandscape  = squeeze(mean(elandscape,1));
meankslandscape = squeeze(mean(kslandscape,2)); 

% 3D energy landscape
figure; 
surf(meanelandscape);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$'); zlabel('Energy');
title(sprintf('%s',modeltype(model)));
xticks(linspace(0,100,5)); xticklabels(-7:3.5:7);
yticks(linspace(0,100,5)); yticklabels(-7:3.5:7);
grid off; shading flat;
set(0,'DefaultAxesColor','none');

% 2D energy landscape
figure;
imagesc(meanelandscape);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
title(sprintf('%s',modeltype(model)));
xticks(linspace(0,100,5)); xticklabels(-7:3.5:7);
xticks([]); yticks([]);
yticks(linspace(0,100,5)); 
yticklabels(-7:3.5:7);
grid off; 
c = colorbar;
c.Label.String = 'Energy'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';

% 2D ks landscapes
ks = string({...
    'Efficiency',...
    'Clustering',...
    'Degree',...
    'Edge length'});

figure;
for k = 1:4
    subplot(2,2,k);
    imagesc(squeeze(meankslandscape(k,:,:)));
    xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
    title(sprintf('%s',ks(k)));
    xticks(linspace(0,100,5)); xticklabels(-7:3.5:7);
    xticks([]); yticks([]);
    yticks(linspace(0,100,5)); 
    yticklabels(-7:3.5:7);
    grid off; 
    c = colorbar;
    c.Label.String = 'KS'; 
    c.Label.Interpreter = 'latex'; 
    c.TickLabelInterpreter = 'latex';
end

%% 7. Plot energy landscapes: all cultures individually

figure; 
set(gcf,'Position',[1000 1000 1250 1500])
for culture = 1:size(MEAenergy,1)
    subplot(4,4,culture);
    sgtitle(sprintf('MEA cultures: %s',modeltype(model)));
    % calculate number of edges
    imagesc(squeeze(elandscape(culture,:,:)));
    n = nnz(squeeze(MEAbinary(:,:,culture)))/2;
    title(sprintf('%g (%g edges)',culture,n));
    xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
    caxis([0 1]);
    xticks([]);
    yticks([]);
end

%% 8. Plot culture wise optimal wiring equations

% find each cultures optimal parameter combination
params_set = [p1(:,model) p2(:,model)];

% visualise
figure;
scatter(params_set(:,1),params_set(:,2),50,'filled','k');
xlim([-7 7]); ylim([-7 7]);
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');

%% 9. Plot specific culture with energy and ks statistics

% set culture
culture = 1;

% make into a graph
g = graph(squeeze(MEAbinary(:,:,culture)));

figure; 
set(gcf,'Position',[1000 1000 1200 300])
n = nnz(squeeze(MEAbinary(:,:,culture)))/2;
sgtitle(sprintf('Culture %g (%g edges)',culture,n));

% plot network
subplot(1,3,1);
h = plot(g,'EdgeAlpha',0.8,'EdgeColor','k','NodeColor','k','Xdata',MEAcoordinates(:,1),'Ydata',MEAcoordinates(:,2),'interpreter','latex');
labelnode(h,1:60,'') % remove node labels
title('Multi-electrode array');
set(gca,'visible','off');

% plot energy landscape
subplot(1,3,2);
imagesc(squeeze(elandscape(culture,:,:)));
title(sprintf('Energy landscape: %s',modeltype(model)));
xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
caxis([0 1]); xticks([]); yticks([]);
c = colorbar;
c.Label.String = 'Energy'; 
c.Label.Interpreter = 'latex'; 
c.TickLabelInterpreter = 'latex';

% plot rule performance
rearrange = [3 2 4 6 7 8 5 9 11 12 13 10 1];
subplot(1,3,3);
h = bar(e(culture,rearrange));
title('Energy distribution of the best network for each generative rule'); 
ylabel('Energy'); 
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]);
xticklabels(modeltype(rearrange));
set(gca,'TickLength',[0 0]);

% 2D ks landscapes
ks = string({...
    'Efficiency',...
    'Clustering',...
    'Degree',...
    'Edge length'});

% plot ks landscape
figure;
sgtitle(sprintf('Culture %g',culture));
for k = 1:4
    subplot(2,2,k);
    imagesc(squeeze(kslandscape(k,culture,:,:)));
    xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
    title(sprintf('%s',ks(k)));
    xticks(linspace(0,100,5)); xticklabels(-7:3.5:7);
    xticks([]); yticks([]);
    yticks(linspace(0,100,5)); 
    yticklabels(-7:3.5:7);
    grid off; 
    caxis([0 1]);
    c = colorbar;
    c.Label.String = 'KS'; 
    c.Label.Interpreter = 'latex'; 
    c.TickLabelInterpreter = 'latex';
end

%% 10. Plot cultures in order of size

figure;
set(gcf,'Position',[1000 1000 2000 300])

% initialise
ind = [length(nnzrank)+1:2*length(nnzrank)];

for n = 1:length(nnzrank);
    
    % take data
    k = nnzrank(n);
    g = graph(squeeze(MEAbinary(:,:,k)));
    
    % plot graph
    subplot(2,length(nnzrank),n);
    h = plot(g,'EdgeAlpha',0.8,'EdgeColor','k','NodeColor','k','Xdata',MEAcoordinates(:,1),'Ydata',MEAcoordinates(:,2),'interpreter','latex');
    labelnode(h,1:60,'');
    xticks([]); yticks([]);
    
    % plot landscape
    subplot(2,length(nnzrank),ind(n));
    imagesc(squeeze(elandscape(k,:,:)));
    xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
    caxis([0 1]); xticks([]); yticks([]);
    
end

%% 11. Plot cultures in order of age (Organoid data only)

% *** set organoid age ***
age = MEAagesCS30;

% sort by age
[agei,agerank] = sort(age);

% initialise
ind = [length(agerank)+1:2*length(agerank)];

for n = 1:length(agerank);
    
    % take data
    k = agerank(n);
    g = graph(squeeze(MEAbinary(:,:,k)));
    
    % plot graph
    subplot(2,length(agerank),n);
    h = plot(g,'EdgeAlpha',0.8,'EdgeColor','k','NodeColor','k','Xdata',MEAcoordinates(:,1),'Ydata',MEAcoordinates(:,2),'interpreter','latex');
    labelnode(h,1:60,'');
    xticks([]); yticks([]);
    
    % plot landscape
    subplot(2,length(agerank),ind(n));
    imagesc(squeeze(elandscape(k,:,:)));
    xlabel('$$ \eta $$'); ylabel('$$ \gamma $$');
    caxis([0 1]); xticks([]); yticks([]);
    
end

% find indices of those with the same age
[c ia ic] = unique(age);

% mean energy and parameters in binned ages
binned_e = [];
binned_p = [];
binned_n = [];
for k = 1:length(c);
    ind = find(ic == k);
    binned_e(k,:) = mean(e(ind,:));
    binned_p(k,:,1) = mean(p1(ind,:));
    binned_p(k,:,2) = mean(p2(ind,:));
    binned_n(k) = length(ind);
end

% plot energy over ages for all rules
figure;
plot(binned_e);
legend(modeltype,'Location','NorthEastOutside');
set(gca,'TickLength',[0 0]);
xticks(1:length(c));
xlabel('Days');
ylabel('Energy');
% form x tick labels
v = {};
for u = 1:length(c);
    v{u} = sprintf('%g n=%g',c(u),binned_n(u));
end
xticklabels(v);

% plot p1 over ages for all rules
figure;
plot(squeeze(binned_p(:,:,1)));
legend(modeltype,'Location','NorthEastOutside');
set(gca,'TickLength',[0 0]);
xticks(1:length(c));
xlabel('Days');
ylabel('$$\eta$$');
% form x tick labels
v = {};
for u = 1:length(c);
    v{u} = sprintf('%g n=%g',c(u),binned_n(u));
end
xticklabels(v);
