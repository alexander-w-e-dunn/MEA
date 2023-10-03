% 
% 
% 
% 
clear all
% cd 'C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data'
% load('3d_matrix_coords_MEAs.mat')
% for i=1:length(coord)
%    channels(i) = str2num([num2str(coord(i,1)),num2str(coord(i,2))]) 
% end
% channels = channels'
% adjM = adjM_all(:,:,1);

cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
% load('MPT200108_1A_DIV14_cSpikes_L-0.0627_RF1_adjM_0.05.mat');
% load('MPT200209_1A_DIV14_cSpikes_L-0.0627_RF1_adjM_0.05.mat');
load('MPT200115_1A_DIV14_cSpikes_L-0.0627_RF1_adjM_0.05.mat');
% load('MPT200115_1B_DIV14_cSpikes_L-0.0627_RF1_adjM_0.05.mat');
% load('MPT200209_1A_DIV14_cSpikes_L-0.0627_RF1_adjM_0.05.mat');
% load('MPT200209_1B_DIV14_cSpikes_L-0.0627_RF1_adjM_0.05.mat');
% load('MPT200209_2A_DIV14_cSpikes_L-0.0627_RF1_adjM_0.05.mat'); %WT
adjM = adjM - eye(size(adjM));
% adjM(adjM>0.5)= 1;
% adjM(adjM~=1) = 0;

adjM = threshold_proportional(adjM, 0.1);
% adjM = threshold_absolute(adjM, 0.5);
adjM = weight_conversion(adjM,'binarize');
adjM = weight_conversion(adjM,'autofix');
adjM(find(channels == 15),:) = 0;
adjM(:,find(channels == 15)) = 0;
% figure
% imagesc(adjM)
% figure
% circularGraph(adjM)
addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'))
MEAgraphics(1)
[C,Q] = community_louvain(adjM,1);     % get community assignments
% [Ci Q] = modularity_und(adjM,1.2);
% [X,Y,INDSORT] = grid_communities(C); % call function
Ci = C;
[X,Y,INDSORT] = grid_communities(Ci); % call function
f2 = figure;
imagesc(adjM(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);             % plot community boundaries
a = gca;
a.FontName = 'Arial';
a.FontSize = 12;
f2.Position = [620   463   552   495];
a.TickLength= 0.01*[1 1]; box off
a.TickDir = 'out';
xticks(1:2:length(adjM)-1);
xtickangle(90);
yticks(2:2:length(adjM));
ytickangle(0);
xticklabels(channels(INDSORT(1:2:length(adjM)-1)));
yticklabels(channels(INDSORT(2:2:length(adjM))));
xlabel('Electrode coordinate (xy)')
ylabel('Electrode coordinate (xy)')
labels = channels(INDSORT);
col_list = lines;
for i = 1 : length(adjM)
    cols(i,1:3) = col_list(Ci(i),:);
end
for i = 1 : length(labels)
    strings{i} = num2str(labels(i));
end

set(groot,'defaultTextFontName','Arial')
set(groot,'defaultTextFontSize',12)

f1 = figure;
circularGraph(adjM(INDSORT,INDSORT),'Colormap',...
    cols(INDSORT,:),'Label',strings)
s = findobj(gcf, 'Type', 'line');
for i =1:length(s)-length(adjM)*2
    s(i).Color = 0.5*[1 1 1];
    %s(i).Color(4) = 0.2;
    s(i).LineWidth = 1;
end

% count = 0;
for i = 1 : length(s)
%     count = count +1 ;
    s(i).MarkerFaceColor = s(i).Color;
%     s(i).Parent.FontName = 'Arial';
%     s(i).Parent.FontSize = 20;
end


