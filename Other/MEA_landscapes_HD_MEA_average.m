
%% load data
load HD_MEA_energy_timepoint_5.mat
%%
modeltype = string({...
    'sptl',...
    'neighbors',...
    'matching',...
    'clu-avg',...
    'clu-min',...
    'clu-max',...
    'clu-diff',...
    'clu-prod',...
    'deg-avg',...
    'deg-min',...
    'deg-max',...
    'deg-diff',...
    'deg-prod'});

f1 = figure;

% for culture = 1 : size(HD_MEA_energy,1)
    
MEAenergy=nanmean(HD_MEA_energy,1);

for rule = 2 : size(MEAenergy,2)
    
    % get lowest energy value for each rule
    %     optimum_params = find
    
    nColumns = 3;
    subplot(ceil(size(MEAenergy,2)/nColumns) , nColumns, rule-1)
    imagesc(reshape(MEAenergy(1,rule,:),sqrt(size(MEAenergy,3)) , sqrt(size(MEAenergy,3)) ))
    set(gca, 'XTick', [1 sqrt(size(MEAenergy,3))], 'XTickLabel', [-7 7])
    set(gca, 'YTick', [1 sqrt(size(MEAenergy,3))], 'YTickLabel', [-7 7])
    ylabel('\eta')
    xlabel('\gamma')
    title(modeltype(rule))
    caxis([0 1])
end
    
% end

f1.Position = [680   137   560   841];

% imagesc(reshape(HDMEA1449energy(1,3,:),10,10));
% xticklabels

%% boxplot

figure
a = min(HD_MEA_energy,[],3);
addpath(genpath('C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts\Other_scripts\RF_scripts'));
addpath(genpath('C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts'));
% notBoxPlot(a)
% notBoxPlotRF(a)
% notBoxPlotRF(a,'jitter',0.5)
% notBoxPlotRF(a,'jitter',1)
% notBoxPlotRF(a,'r',2)
% notBoxPlotRF(a,[0.5 0.5 0.5],2)
% 
% notBoxPlotRF(a,{[0 0 0], cols(1), cols(1), cols(2), cols(2), cols(2), cols(2), cols(2), cols(3), cols(3), cols(3), cols(3), cols(3)},2)
% length([1 1 1; 1 1 1]
% length([1 1 1; 1 1 1])
% [1 1 1; 1 1 1]
% ones(3,13)
% notBoxPlotRF(a,ones(3,13),2)
% notBoxPlotRF(a,0.5*ones(3,13),2)
% notBoxPlotRF(a,1,2)
% notBoxPlotRF(a,1:13,1,2)
% notBoxPlotRF(a,1:13,[0 0 0],2)
% notBoxPlotRF(a,1:13,[0 0 0],1)
% hold off
% notBoxPlotRF(a,1:13,[0 0 0],1)
% notBoxPlotRF(a,1:13,[0 0 0],0.5)
% figure
% notBoxPlotRF(a,1:13,[0 0 0],0.5)
% figure
% notBoxPlotRF(a,1:13,[0 0 0],1)
% notBoxPlotRF(a,1:13,[0 0 0; 1 0 0],1)
% figure
% notBoxPlotRF(a,1,[0 0 0],1)
% notBoxPlotRF(a(:,1),1,[0 0 0],1)
% hold on
% notBoxPlotRF(a(:,2),2,[0 0 0],1)
% notBoxPlotRF(a(:,[2 3]),[2 3],[0 0 0],1)
% notBoxPlotRF(a(:,[2 3]),[2 3],'r',1)
% notBoxPlotRF(a(:,[2 3]),[2 3],cols(1),1)
% cols(1)
% cols(:,1)
% notBoxPlotRF(a(:,[2 3]),[2 3],cols(1,:),1)
% figure
% notBoxPlotRF(a,1,[0 0 0],1)
% notBoxPlotRF(a(:,1),1,[0 0 0],1)
% notBoxPlotRF(a(:,[2 3]),[2 3],cols(1,:),1)
% notBoxPlotRF(a(:,[2 3]),[4:8],cols(2,:),1)
% notBoxPlotRF(a(:,[4:8]),[4:8],cols(2,:),1)
% hold on
% notBoxPlotRF(a(:,[9:13]),[9:13],cols(3,:),1)
% xlim([0 14])
% figure
% notBoxPlotRF(a,1,[0 0 0],1,'style','sdline')
% notBoxPlotRF(a(:,1),1,[0 0 0],1,'style','sdline')
% figure
% notBoxPlotRF(a(:,1),1,[0 0 0],1,'style','line')
% figure
% notBoxPlotRF(a(:,1),1,[0 0 0],1,'style','sdline')
% notBoxPlotRF(a(:,[2 3]),[2 3],cols(1,:),1,'style','sdline')

cols=lines;
figure
notBoxPlotRF(a(:,1),1,[0 0 0],1,'style','sdline')
hold on
notBoxPlotRF(a(:,[2 3]),[2 3],cols(1,:),1,'style','sdline')
hold on
notBoxPlotRF(a(:,[4:8]),[4:8],cols(5,:),1,'style','sdline')
hold on
notBoxPlotRF(a(:,[9:13]),[9:13],cols(3,:),1,'style','sdline')
xlabel('Rule')
xlim([0 14])
xticks(1:13)
ylabel('Energy')
b=gca;
b.TickDir = 'out';
b.FontSize = 16;
aesthetics

xticklabels(modeltype)
xtickangle(45)