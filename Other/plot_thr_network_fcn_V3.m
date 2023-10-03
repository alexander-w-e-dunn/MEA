% script to plot the network and its mst
% Written by AWE Dunn Cambridge, 2020
% note, to use this with MEA data you need to set up the coordinates. this
% involves reversing the y coordinate
% 
% INPUTS:   fileName = .mat file name containing the matrix called adjM and a
%                       set of xy coordinates saved in a vector variable called channels such
%                       that channels(i) are the cordinates xy (currently only works for x = 0-9
%                       and y = 0 - 9
%           thrOption = 'proportional' or 'absolute'
%           thr = v     alue for thresholding as a proportion of connections
%                       to keep or an absolute value below which edges are removed
%           type  =     'weighted' or 'binary'
%           xy =        matrix in which the first column is a vector of x
%                       coordinates and the second column is a vector of y coords
% 
% OUTPUTS: will save PNG to working directory of the network and then one
% 
% example use:
% 
% fileName = '200127_FTDOrg_GrpD_5B_Slice7_mSpikes_3_adjM_0.175.mat';
% thrOption = 'proportional';
% thr = 0.4;
% type = 'weighted'
% plot_thr_network_fcn(fileName,thrOption,thr,type)

% clear all
% set(groot,'defaultAxesFontName','Arial')
% set(groot,'defaultAxesFontSize',22)
% set(groot,'defaultAxesLineWidth',2)
% set(groot,'defaultFigurePosition',[700   435   520   350])

function plot_thr_network_fcn_V3(adjM,thrOption,thr,type,xy)
%% get the adjM, threshold, binarise and plot
MEAgraphics(1);

% take top 5 connections
KadjM = weight_conversion(adjM, 'autofix');
for cell = 1 : length(adjM)
    Ksorted = sort(KadjM(:,cell),'descend');
    klevel = Ksorted(5);
    vec_to_remove = find(KadjM(:,cell) < klevel);
    adjM(vec_to_remove,cell) = 0;
end

% load(fileName,'adjM','channels');
% remove ref
% adjM(find(channels == 15),:) = 0;
% adjM(:,find(channels == 15)) = 0;
% remove negative connections for sttc
adjM(adjM < 0) = 0; 
% take absolute value to include negative correlations if not using sttc
% adjM = abs(adjM);
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
xc = round(xy(:,1),0);
yc = round(xy(:,2),0);
% get coords for plotting
[X,Y] = adjacency_plot_und(badjM,[xc, yc]);
% X = weight_conversion(X, 'autofix'); Y = weight_conversion(Y, 'autofix');
f1 = figure; f1.Position = [700   250   792   535];
% plot edges
plot(X,Y,'Visible','off'); aesthetics; axis off;hold on
max_ew = 3; % maximum edge width for plotting
min_ew = 0.5;     % min edge width

if strcmp(type,'weighted')
    if strcmp(thrOption,'absolute')
        Emin = thr;
    else
        Emin = min(badjM1(badjM == 1));
    end
end

maxX = max(xc) + 1;
maxY = max(yc) + 1;

if strcmp(type,'weighted')
    for elecA = 1:length(xy)
        for elecB = 1:length(xy)
            if badjM1(elecA,elecB) >= Emin & elecA ~= elecB;
                p = plot([xc(elecA),xc(elecB)],[yc(elecA),yc(elecB)],'LineWidth',...
                    min_ew + (max_ew-min_ew)*((badjM1(elecA,elecB)-Emin)/(1-Emin)),'Color',0.5*[1 1 1]);
                % set transparency of lines
                p.Color(4) = 0.2;
            end
        end
    end
elseif strcmp(type,'binary')
    for elecA = 1:length(xy)
        for elecB = 1:length(xy)
            if badjM1(elecA,elecB) >= Emin & elecA ~= elecB;
                p = plot([xc(elecA),xc(elecB)],[yc(elecA),yc(elecB)],'LineWidth',...
                    1,'Color',0.5*[1 1 1]);
                % set transparency of lines
                p.Color(4) = 0.2;
            end
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
    plot(xc(i),yc(i),'ok','MarkerSize',5,...
        'MarkerEdgeColor',colours(i,:),...
        'LineWidth',1)
end

% add hidden legend data and adjust offset 
nodemin = 2; nodemid = 6; nodemax = 11; 
L1 = plot([  10000],[  10000],'ok','MarkerSize',5*sqrt([1 + nodemin ]),...
    'MarkerFaceColor','k');
L2 = plot([  10000],[  10000],'ok','MarkerSize',5*sqrt([1 + nodemid ]),...
    'MarkerFaceColor','k');
L3 = plot([  10000],[  10000],'ok','MarkerSize',5*sqrt([1 + nodemax ]),...
    'MarkerFaceColor','k');

xlim([0 maxX+1]); ylim([0 maxY]);

if strcmp(type,'binary')
    L = legend([ L1 L2 L3], {num2str(nodemin-1),num2str(nodemid-1),num2str(nodemax-1)},'box','off');
    L.Position = [0.8800    0.45   0.1364    0.3598];
end

% colorbar for edgeweight
if strcmp(type,'weighted')
    cb1 = colorbar;
    cb1.Position = [ 0.93    0.1095    0.020    0.15 ];
    cb1.Ticks = [0 1];
end
hold on

% add scale bar
plot([0:200],zeros(201,1),'-k')
text(0,-maxY*0.02,'200 \muM','FontSize',12);

% add hidden legend data for edge weight
if strcmp(type,'weighted')
    if strcmp(thrOption,'absolute')
        Emin = thr;
    else
        Emin = min(badjM1(badjM == 1));
    end
    % added correction here if minimum edge weight = 1
    if Emin == 1 
        Emin = Emin-0.0001;
    elseif isempty(Emin)
        Emin = 0;
    end
    Emed = Emin+((1-Emin)/2); %middle of min and max
    Emax = 1; %max value an edge weight can be
    l4 = plot([1 2],[-10 -10],'LineWidth',...
        min_ew + (max_ew-min_ew)*((Emin -Emin)/(1-Emin)),'Color',0.5*[1 1 1]);
    l5 = plot([1 2],[-10 -10],'LineWidth',...
        min_ew + (max_ew-min_ew)*((Emed -Emin)/(1-Emin)),'Color',0.5*[1 1 1]);
    l6 = plot([1 2],[-10 -10],'LineWidth',...
        min_ew + (max_ew-min_ew)*((Emax -Emin)/(1-Emin)),'Color',0.5*[1 1 1]);
    LX1 = plot(1,-1,'o','MarkerSize',0.000000001,'MarkerEdgeColor','w','MarkerFaceColor','w');
    LX2 = plot(1,-1,'o','MarkerSize',0.000000001,'MarkerEdgeColor','w','MarkerFaceColor','w');
    ylim([0 maxY])
    % add legend
%     L = legend([ L1 L2 L3 LX1 LX2 l4 l5 l6], ...
%         {num2str(nodemin-1),num2str(nodemid-1),num2str(nodemax-1),[''],[''],...
%         num2str(Emin,'%.2f') , num2str(Emed,'%.2f') , num2str(Emax,'%.2f')},'box','off');
    L = legend([ LX1 LX1 LX1 LX1 LX2 l4 l5 l6], ...
        {num2str(nodemin-1),num2str(nodemid-1),num2str(nodemax-1),[''],[''],...
        num2str(Emin,'%.2f') , num2str(Emed,'%.2f') , num2str(Emax,'%.2f')},'box','off');

    L.Position = [0.8800    0.45   0.1364    0.3598];
end

% add titles of legends

% text(8.5,7.8,'degree','FontSize',22);
% if strcmp(type,'weighted')
%     text(8.5,5,'weight','FontSize',12);
%     text(8.5,2.2,'strength','FontSize',12);
% end
f1.Position = [871   177   742   608];
% saveas(f1,strcat(fileName,'_network_plot_',type,'_',thrOption,...
%     '_',num2str(thr),'.png'))
% close(f1)
MEAgraphics(0);
end
