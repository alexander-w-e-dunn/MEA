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
% 
% OUTPUTS: will save PNG to working directory of the process of
% thresholding and binarising
% 
% example use:
% 
% fileName = '200127_FTDOrg_GrpD_5B_Slice7_mSpikes_3_adjM_0.175.mat';
% thrOption = 'proportional';
% thr = 0.4;
% plot_binary_approach(fileName,thrOption,thr)

% clear all
% set(groot,'defaultAxesFontName','Arial')
% set(groot,'defaultAxesFontSize',22)
% set(groot,'defaultAxesLineWidth',2)
% set(groot,'defaultFigurePosition',[700   435   520   350])

function plot_binary_approach(fileName,thrOption,thr)
%% get the adjM, threshold, binarise and plot
load(fileName,'adjM','channels');
% remove negative connections for sttc
adjM(adjM < 0) = 0;

f2=figure;
f2.Position = [  614         271        1240         708  ];

% plot  full adjM
subplot(2,3,1);imagesc(adjM); xticklabels('');yticklabels('');caxis([0 1]);
ylabel('electrode');xlabel('electrode');

% take absolute value to include negative correlations if not using sttc
% adjM = abs(adjM);
% threshold the matrix
if strcmp(thrOption,'proportional')
    badjM1 = threshold_proportional(adjM, thr);
elseif strcmp(thrOption,'absolute')
    badjM1 = threshold_absolute(adjM, thr);
end
badjM1 = weight_conversion(badjM1, 'autofix');
% plot  thresholded wweighted adjM
subplot(2,3,2);imagesc(badjM1); xticklabels('');yticklabels('');
ylabel('');xlabel('electrode');caxis([0 1]);
% get sum of weights for each node > threshold
w = sum(badjM1)'; c = colormap;
% binarise matrix
badjM = weight_conversion(badjM1, 'binarize');
% plot  thresholded wweighted adjM
subplot(2,3,3);imagesc(badjM); xticklabels('');yticklabels('');
ylabel('');xlabel('electrode');caxis([0 1]);
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
% plot edges
types = {'weighted','binary'};
plotn = [5 6];
for i = 1:length(types)
    subplot(2,3,plotn(i))
    type = types{i};
    plot(X,Y,'Visible','off'); aesthetics; axis off;hold on
    max_ew = 3; % maximum edge width for plotting
    min_ew = 0.5;     % min edge width
    
    if strcmp(type,'weighted')
        for elecA = 1:length(channels)
            for elecB = 1:length(channels)
                if badjM1(elecA,elecB) >= thr & elecA ~= elecB;
                    p = plot([xc(elecA),xc(elecB)],[9-yc(elecA),9-yc(elecB)],'LineWidth',...
                        min_ew + (max_ew-min_ew)*((badjM1(elecA,elecB)-thr)/(1-thr)),'Color',0.5*[1 1 1]);
                    % set transparency of lines
                    p.Color(4) = 0.2;
                end
            end
        end
    elseif strcmp(type,'binary')
        for elecA = 1:length(channels)
            for elecB = 1:length(channels)
                if badjM1(elecA,elecB) >= thr & elecA ~= elecB;
                    p = plot([xc(elecA),xc(elecB)],[9-yc(elecA),9-yc(elecB)],'LineWidth',...
                        max_ew,'Color',0.5*[1 1 1]);
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
        plot(xc(i),9-yc(i),'ok','MarkerSize',2*sqrt(1+degree(i)),...
            'MarkerFaceColor',colours(i,:))
    end
        
end

%%colorbar for edgeweight
subplot(2,3,4)
% add hidden legend data and adjust offset 
nodemin = 2; nodemid = 6; nodemax = 11; 
L1 = plot([  10],[  10],'ok','MarkerSize',5*sqrt([1 + nodemin ]),...
    'MarkerFaceColor','k'); hold on
L2 = plot([  10],[  10],'ok','MarkerSize',5*sqrt([1 + nodemid ]),...
    'MarkerFaceColor','k');
L3 = plot([  10],[  10],'ok','MarkerSize',5*sqrt([1 + nodemax ]),...
    'MarkerFaceColor','k');

xlim([1 9]); ylim([0 8]);

% add hidden legend data for edge weight
Emin = min(badjM1(badjM == 1));
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
xlim([1 9]); ylim([0 8]);

% add legend
L = legend([ L1 L2 L3 LX1 l4 l5 l6], ...
    {num2str(nodemin-1),num2str(nodemid-1),num2str(nodemax-1),[''],...
    num2str(Emin,'%.2f') , num2str(Emed,'%.2f') , num2str(Emax,'%.2f')},'box','off');
axis off
L.Position = [0.2519    0.1191    0.0831    0.3164];

cb1 = colorbar;
cb1.Location = 'westoutside';
cb1.Label.String = 'strength';
cb1.Label.FontSize = 22;

% add titles of legends
text(4.7,8.1,'degree','FontSize',22);
text(4.7,4,'weight','FontSize',22);



saveas(f2,strcat(fileName,'_approach_',thrOption,...
    '_',num2str(thr),'.png'))
close(f2)


end
