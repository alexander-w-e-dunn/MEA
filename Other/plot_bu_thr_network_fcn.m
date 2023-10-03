% script to plot the network and its mst
% Written by AWE Dunn Cambridge, 2020
% 
% INPUTS:   fileName = .mat file name containing the matrix called adjM and a
%                       set of xy coordinates saved in a vector variable called channels such
%                       that channels(i) are the cordinates xy (currently only works for x = 0-9
%                       and y = 0 - 9
%           thrOption = 'proportional' or 'absolute'
%           thr = v     alue for thresholding as a proportion of connections
%                       to keep or an absolute value below which edges are removed
%           type  =     'weighted' or 'binary'
% 
% OUTPUTS: will save PNG to working directory of the network and then one
% 
% example use:
% 
fileName = '200127_FTDOrg_GrpD_5B_Slice7_mSpikes_3_adjM_0.175.mat';
thrOption = 'proportional';
thr = 0.4;
plot_wu_thr_network_fcn(fileName,thrOption,thr,type)
type = 'binary'

% clear all
% set(groot,'defaultAxesFontName','Arial')
% set(groot,'defaultAxesFontSize',22)
% set(groot,'defaultAxesLineWidth',2)
% set(groot,'defaultFigurePosition',[700   435   520   350])

function plot_wu_thr_network_fcn(fileName,thrOption,thr,type)
%% get the adjM, threshold, binarise and plot
load(fileName,'adjM','channels');
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
plot(X,Y,'Visible','off'); aesthetics; axis off;hold on
max_ew = 2; % maximum edge width for plotting
min_ew = 0.5;     % min edge width
for elecA = 1:length(channels)
    for elecB = 1:length(channels)
        if badjM1(elecA,elecB) >= thr & elecA ~= elecB;
            p = plot([xc(elecA),xc(elecB)],[yc(elecA),yc(elecB)],'LineWidth',...
                min_ew + (max_ew-min_ew)*((badjM1(elecA,elecB)-thr)/(1-thr)),'Color',0.5*[1 1 1]);
            % set transparency of lines
            p.Color(4) = 0.2; 
        end
    end
end

% convert edges with no weight to slightly above 0
% this is so that their colour can be colour 1 of the colormap
w(w<=0)=0.00001;
% make colormap minimum black rather than blue so that disconnected nodes
% are in black
c(1,:) = 0;
% plot nodes
if strcmp(type,'weighted')
    colours = c(ceil(w*length(c)),:);
else 
    colours = zeros(length(adjM),3)  + 0;
end
for i = 1:length(xc)
    plot(xc(i),yc(i),'ok','MarkerSize',5*sqrt(1+degree(i)),...
        'MarkerFaceColor',colours(i,:))
end

% add hidden legend data and adjust offset 
nodemin = 1; nodemid = 10; nodemax = 50; 
L1 = plot([  10],[  10],'ok','MarkerSize',2*sqrt([1 + nodemin ]),...
    'MarkerFaceColor','k');
L2 = plot([  10],[  10],'ok','MarkerSize',2*sqrt([1 + nodemid ]),...
    'MarkerFaceColor','k');
L3 = plot([  10],[  10],'ok','MarkerSize',2*sqrt([1 + nodemax ]),...
    'MarkerFaceColor','k');

xlim([1 9]);

if strcmp(type,'binary')
    L = legend([ L1 L2 L3], {num2str(nodemin),num2str(nodemid),num2str(nodemax)},'box','off');
    L.Position = [0.8600    0.45   0.1364    0.3598];
end

% colorbar for edgeweight
if strcmp(type,'weighted')
    cb1 = colorbar;
    cb1.Position = [ 0.9    0.1095    0.020    0.15 ];
    cb1.Ticks = [0 1];
end

% add hidden legend data for edge weight
if strcmp(type,'weighted')
    Emin = thr; %min value an edge can be (for legend)
    Emed = thr+((1-thr)/2); %middle of min and max
    Emax = 1; %max value an edge weight can be
    l4 = plot([1 2],[-10 -10],'LineWidth',...
        min_ew + (max_ew-min_ew)*((Emin -thr)/(1-thr)),'Color',0.5*[1 1 1]);
    l5 = plot([1 2],[-10 -10],'LineWidth',...
        min_ew + (max_ew-min_ew)*((Emed -thr)/(1-thr)),'Color',0.5*[1 1 1]);
    l6 = plot([1 2],[-10 -10],'LineWidth',...
        min_ew + (max_ew-min_ew)*((Emax -thr)/(1-thr)),'Color',0.5*[1 1 1]);
    LX1 = plot(1,-1,'o','MarkerSize',0.000000001,'MarkerEdgeColor','w','MarkerFaceColor','w');
    LX2 = plot(1,-1,'o','MarkerSize',0.000000001,'MarkerEdgeColor','w','MarkerFaceColor','w');
    ylim([0 8])
    % add legend
    L = legend([ L1 L2 L3 LX1 LX2 l4 l5 l6], ...
        {num2str(nodemin),num2str(nodemid),num2str(nodemax),[''],[''],...
        num2str(Emin,'%.2f') , num2str(Emed,'%.2f') , num2str(Emax,'%.2f')},'box','off');
    L.Position = [0.8600    0.45   0.1364    0.3598];
end
% add titles of legends
text(8.5,7.8,'degree','FontSize',22);
if strcmp(type,'weighted')
    text(8.5,5,'weight','FontSize',22);
    text(8.5,2.2,'strength','FontSize',22);
end
f1.Position = [871   177   742   608];
saveas(f1,strcat(fileName,'_network_plot_nodalStrength.png'))
close(f1)

end
