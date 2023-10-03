% script to plot the network and its mst
% Written by AWE Dunn and D Akarca, Cambridge, 2020
% 
% INPUTS: fileName: .mat file name containing the matrix called adjM and
% list of xy coordinates as a vector e.g. [12 13 24] for nodes 1 to 3 mean
% node 1 x = 1 and y = 2
%         colour scheme: this is 'colour' for coloured 'black' for
%         greyscale
% 
% OUTPUTS: will save PNG to working directory of the network and then one
% of this compared to its min span tree
% 
% 
% 
% 
% 
clear all
set(groot,'defaultAxesFontName','Arial')
set(groot,'defaultAxesFontSize',22)
set(groot,'defaultAxesLineWidth',2)
set(groot,'defaultFigurePosition',[700   435   520   350])

% file = '200127_FTDOrg_GrpD_5B_Slice7_mSpikes_3_adjM_0.175.mat';
% colourscheme = 'colour';
% thrOption = 'proportional';
% thr = 0.4;

function plot_mst(fileName,colourscheme,thrOption,thr)
%% get the adjM, threshold, binarise and plot
load(file,'adjM','channels');
% threshold the matrix
if strcmp(thrOption,'proportional')
    badjM1 = threshold_proportional(adjM, thr);
elseif strcmp(thrOption,'absolute')
    badjM1 = threshold_absolute(adjM, thr);
end
badjM1 = weight_conversion(badjM1, 'autofix');
% get sum of weights for each node > threshold
w = sum(badjM1)'; c = colormap;
% binarise matrix
badjM = weight_conversion(badjM1, 'binarize');
% calculate degree of each node
degree = sum(badjM)';
% convert sum of weights to mean weight including only edges > threshold
w = w./degree; w(isnan(w)) = 0;
% get node coords
xc = floor(channels/10);
yc = channels - xc*10;
% get coords for plotting
[X,Y] = adjacency_plot_und(badjM,[xc, yc]);
% X = weight_conversion(X, 'autofix'); Y = weight_conversion(Y, 'autofix');
f1 = figure; f1.Position = [700   250   792   535];
% plot edges
p = plot(X,Y,'k'); aesthetics; axis off;hold on
% set transparency of lines
p.Color(4) = 0.2; 
% convert edges with no weight to slightly above 0
% this is so that their colour can be colour 1 of the colormap
w(w<=0)=0.00001;
% make colormap minimum black rather than blue so that disconnected nodes
% are in black
c(1,:) = 0;
% plot nodes
if strcmp(colourscheme,'colour')
    colours = c(ceil(w*length(c)),:);
else 
    colours = zeros(length(adjM),3)  + 0;
end
for i = 1:length(xc)
    plot(xc(i),yc(i),'ok','MarkerSize',5*sqrt(1+degree(i)),...
        'MarkerFaceColor',colours(i,:))
end

% add legend and adjust offset 
nodemin = 1; nodemid = 10; nodemax = 50; 
L1 = plot([  10],[  10],'ok','MarkerSize',2*sqrt([1 + nodemin ]),...
    'MarkerFaceColor','k');
L2 = plot([  10],[  10],'ok','MarkerSize',2*sqrt([1 + nodemid ]),...
    'MarkerFaceColor','k');
L3 = plot([  10],[  10],'ok','MarkerSize',2*sqrt([1 + nodemax ]),...
    'MarkerFaceColor','k');

L = legend([ L1 L2 L3], {num2str(nodemin),num2str(nodemid),num2str(nodemax)},'box','off');
L.Position = [0.86    0.32    0.1364    0.1440];

xlim([0 9]);

if strcmp(colourscheme,'colour')
    cb1 = colorbar;
    cb1.Position = [ 0.9    0.1095    0.020    0.15 ];
    cb1.Ticks = [0 1];
end
saveas(f1,strcat(fileName,'_netplot.png'))
close(f1)

%% plot adjM and its min span tree together
% for i = 1:10
f2=figure;
f2.Position = [  680   558   572   420];

% plot binary full adjM
subplot(2,2,1);imagesc(badjM); xticklabels('');yticklabels('');

% plot full network
subplot(2,2,3)
p = plot(X,Y,'k'); aesthetics; axis off;hold on
p.Color(4) = 0.2;

if strcmp(colourscheme,'colour')
    colours = c(ceil(w*length(c)),:);
else 
    colours = zeros(length(adjM),3)  + 0;
end
for i = 1:length(xc)
    plot(xc(i),yc(i),'ok','MarkerSize',2*sqrt(1+degree(i)),...
        'MarkerFaceColor',colours(i,:))
end
xlim([0 9]);
ylim([0 8])

% % get weighted back bone 
% [CIJtreew] = mst_wu(badjM1);
% % get weights
% w2 = sum(CIJtreew)';
% getbinary backbone
[CIJtree] = mst_wu(badjM);
% get weights
w_vec = find(CIJtree == 1);
CIJtree_w = zeros(size(CIJtree)); 
CIJtree_w(w_vec) = badjM1(w_vec);
w2 = sum(CIJtree_w)';
% get metrics
degree2 = sum(CIJtree)';
w2 = w2./degree2; w2(isnan(w2)) = 0;
% remove 0s from weights for the purpose of plotting (0.0001 will be
% rounded up to 1 so that the colour of the node will be row 1 of the
% colormap
w2(w2 == 0) = 0.00001;
% plot adjM of backbone
subplot(2,2,2);imagesc(CIJtree);xticklabels('');yticklabels('');
% plot backbone network
subplot(2,2,4)
% get coords for mst
[X,Y] = adjacency_plot_und(CIJtree,[xc, yc]);
% plot edges
p = plot(X,Y,'k'); aesthetics; axis off;hold on
% set transparency of lines
p.Color(4) = 0.2;
% plot nodes
if strcmp(colourscheme,'colour')
    colours = c(ceil(w2*length(c)),:);
else 
    colours = zeros(length(adjM),3)  + 0;
end

for i = 1:length(xc)
    plot(xc(i),yc(i),'ok','MarkerSize',2*sqrt(1+degree2(i)),...
        'MarkerFaceColor',colours(i,:))
end
% add legend and adjust offset 
nodemin = 1; nodemid = 10; nodemax = 50; 
L1 = plot([  10],[  10],'ok','MarkerSize',2*sqrt([1 + nodemin ]),...
    'MarkerFaceColor','k');
L2 = plot([  10],[  10],'ok','MarkerSize',2*sqrt([1 + nodemid ]),...
    'MarkerFaceColor','k');
L3 = plot([  10],[  10],'ok','MarkerSize',2*sqrt([1 + nodemax ]),...
    'MarkerFaceColor','k');

L = legend([ L1 L2 L3], {num2str(nodemin),num2str(nodemid),num2str(nodemax)},'box','off');
L.Position = [0.86    0.32    0.1364    0.1440];

xlim([0 9]);

% add colorbar
if strcmp(colourscheme,'colour')
    cb = colorbar;
    cb.Position = [ 0.9    0.1095    0.020    0.15 ];
    cb.Ticks = [0 1];
end

saveas(f2,strcat(fileName,'_netplotVSminSpanTree.png'))
close(f2)

end
