% script to plot in MEA layout the graph metrics

%{

INPUTS:
z - the metric to colour e.g. edge weight 
z2 - the metric to size the circles. e.g. node degree

future updates:
(need to make a proportion of colour thn relabel
colorbar (e.g. mean connection weight) 

%}
close all; clear all
% simulate data
cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
load('MPT190403_6B_DIV28_cSpikes_L0.mat','channels')
coords = num2str(channels);
x = str2num(coords(1:end,1));
y = 9 - str2num(coords(1:end,2));

%need to add section that calibrates the z values (node sizes) to values
%between 1 and 25
%need to add section that calbriates the z2 values to values between 0 and
%1; then I need to plot and then do the inverse of this calibration on the
%colorbar and legend e.g. edit cb.TickLabels

z = round(abs(rand(60,1)*25));
% z(find(z==0))=1;
z(find(channels == 15)) = 0;

z2 = abs(rand(60,1)); 
z2(find(z2 < 0)) = 0;

%option 1: create colour vector that colours red or blue depending if in
%rich club; requires adding a legend
%option 2: colour according to another metric; requires adding a colorbar

%get colormap

F1 = figure;
F1.OuterPosition = [750   533   456   467];
hold on

mycolours = colormap;
legdata = ['01'; '05'; '10'];
l1 = plot(1,-1,'o','MarkerSize',str2num(legdata(1,:)),'MarkerEdgeColor','k');
l2 = plot(2,-1,'o','MarkerSize',str2num(legdata(2,:)),'MarkerEdgeColor','k');
l3 = plot(3,-1,'o','MarkerSize',str2num(legdata(3,:)),'MarkerEdgeColor','k');

for i = 1:length(channels)
    if z(i)>0
    plot(x(i),y(i),'o','MarkerSize',z(i),'MarkerFaceColor',...
        mycolours(ceil(length(mycolours)*z2(i)),1:3),'LineStyle','none',...
        'MarkerEdgeColor','k') ;
    % where x and y are MEA coords and z is graph metric value
    end
end

set(gca,'color','none')
ylim([0 9])
xlim([0 9])

a = gca;
xticks(1:8)
yticks(1:8)
yvalues = a.YTickLabels ;
yvalues = sort(str2num(cell2mat(yvalues)),'descend');
a.YTickLabels = yvalues;


% % set(a,'xcolor','none')
% a.XAxis.Label.Color=[0 0 0];
% a.XAxis.Label.Visible='on';
%
set(a, 'XAxisLocation', 'top')
a.XAxis.FontName = 'Arial';
a.XAxis.FontSize = 14;
a.YAxis.FontName = 'Arial';
a.YAxis.FontSize = 14;
% a.XAxis.Color = [0.5 0.7 0.8]
% a.YAxis.Color = [0.5 0.7 0.8]

cb = colorbar;
cb.FontName = 'Arial';
cb.FontSize = 14;
cb.FontName = 'Arial';
cb.FontSize = 14;
cb.Label.String = 'Mean edge strength';
cb.Label.FontSize = 14;
cb.Label.FontName = 'Arial';
cb.Label.FontWeight = 'bold';
cb.Label.FontSize = 14;

% need to add colorbar labelling for node color

%may require having colorbar at the bottom
% cb.Location = 'Southoutside'

% need to add legend for node sizes

l4 = legend([l1 l2 l3],'Orientation','horizontal','Position',...
    [0.3 0.07 0 0],'Box','off','String',legdata,...
    'FontSize',14);
l4.Title.String = 'Node degree';
l4.Title.FontSize = 14;
l4.Title.FontName = 'Arial';
l4.Title.FontWeight = 'bold';

for j=1:length(cb.TickLabels)
    cbtlabels(j) = str2num(cb.TickLabels{j});
end






