
% script for comparing groups

%get data
    % lab PC
% cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
    % hp laptop
cd 'C:\Users\owner\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Data'
load('manuel_3_organoidDev_April15.mat')
%back to scripts directory
    % lab PC
% cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
    % hp laptop
addpath(genpath('C:\Users\owner\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts'))

MEAgraphics(1);disp('graphics set')
%% betweenness centrality
clear groups
for i = 1:length(output)
    groups(i) = output(i).grp;
end
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).BC_n;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).BC_n;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).BC_n;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).BC_n;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).BC_n;
end

close all; notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('normalised betweenness centraility')

%% degree
figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).meanDegree;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).meanDegree;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).meanDegree;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).meanDegree;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).meanDegree;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('average node degree')


%% nrmlsd_cnnctvty

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).nrmlsd_cnnctvty;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).nrmlsd_cnnctvty;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).nrmlsd_cnnctvty;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).nrmlsd_cnnctvty;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).nrmlsd_cnnctvty;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('normalised average functional connectivity')

%% meanSTTC

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).meanSTTC;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).meanSTTC;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).meanSTTC;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).meanSTTC;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).meanSTTC;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('average functional connectivity')

%% firing rate


figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).mean_FR;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).mean_FR;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).mean_FR;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).mean_FR;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).mean_FR;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('average firing rate (spikes/s)')

%% netw_density

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).netw_density;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).netw_density;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).netw_density;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).netw_density;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).netw_density;
end

notBoxPlot([dataE*100],5,'style','sdline')
hold on
notBoxPlot([dataD*100],4,'style','sdline')
notBoxPlot([dataC*100],3,'style','sdline')
notBoxPlot([dataB*100],2,'style','sdline')
notBoxPlot([dataA*100],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Network density (%)')

%% netw_size

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).netw_size;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).netw_size;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).netw_size;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).netw_size;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).netw_size;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Network size (# nodes)')

%% CC


figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).CC;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).CC;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).CC;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).CC;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).CC;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Normalised clustering')

%% PL

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).PL;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).PL;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).PL;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).PL;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).PL;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Normalised characteristic path length')
%% SW
clear groups
for i = 1:length(output)
    groups(i) = output(i).grp;
end
ind = find(groups =='E');

F1 = figure;
F1.OuterPosition = [750   533   456   467];

ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).SW2;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).SW2;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).SW2;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).SW2;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).SW2;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Normalised small world coefficient')
ylim([0.8 2.5])
yticks([1 1.5 2 2.5])
hold on
plot([0:6],ones(1,7),'LineStyle',':','Color','k','LineWidth',2)
a = gca;
a.XAxis.FontName = 'Arial';
a.XAxis.FontSize = 14;
a.YAxis.FontName = 'Arial';
a.YAxis.FontSize = 14;

%% RC

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).RC;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).RC;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).RC;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).RC;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).RC;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Rich club coefficient')

%% # RC nodes

clear groups
for i = 1:length(output)
    groups(i) = output(i).grp;
end
ind = find(groups =='E');

F1 = figure;
F1.OuterPosition = [750   533   456   467];

ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).num_RC_nodes_mean;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Count of rich club nodes')
ylim([0 14])
% yticks([1 1.5 2 2.5])
hold on
a = gca;
a.XAxis.FontName = 'Arial';
a.XAxis.FontSize = 14;
a.YAxis.FontName = 'Arial';
a.YAxis.FontSize = 14;

%%  +AGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FR over age 

% output2=output;
% [~,index] = sortrows([output2.mean_FR].'); output2 = output2(index); clear index
% output2(1) = [];
% output2(78:79) = [];
clearvars -except output output2

for i = 1:length(output2)
    groups(i) = output2(i).age_DIV;
end
groups1 = unique(groups);

F1 = figure;
F1.OuterPosition = [750   533   456   467];

ind = find(groups ==groups1(1));

for i = 1:length(ind)
    data1(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(2));
for i = 1:length(ind)
    data2(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(3));
for i = 1:length(ind)
    data3(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(4));
for i = 1:length(ind)
    data4(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(5));
for i = 1:length(ind)
    data5(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(6));
for i = 1:length(ind)
    data6(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(7));
for i = 1:length(ind)
    data7(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(8));
for i = 1:length(ind)
    data8(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(9));
for i = 1:length(ind)
    data9(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(10));
for i = 1:length(ind)
    data10(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(11));
for i = 1:length(ind)
    data11(i) = output2(ind(i)).mean_FR;
end

clear ind; ind = find(groups ==groups1(12));
for i = 1:length(ind)
    data12(i) = output2(ind(i)).mean_FR;
end

notBoxPlot([data1],groups1(1),'style','sdline')
hold on
notBoxPlot([data2],groups1(2),'style','sdline')
notBoxPlot([data3],groups1(3),'style','sdline')
notBoxPlot([data4],groups1(4),'style','sdline')
notBoxPlot([data5],groups1(5),'style','sdline')
notBoxPlot([data6],groups1(6),'style','sdline')
notBoxPlot([data7],groups1(7),'style','sdline')
notBoxPlot([data8],groups1(8),'style','sdline')
notBoxPlot([data9],groups1(9),'style','sdline')
notBoxPlot([data10],groups1(10),'style','sdline')
notBoxPlot([data11],groups1(11),'style','sdline')
notBoxPlot([data12],groups1(12),'style','sdline')
% ylim([10 60])
aesthetics

%% # RC nodes over age 

% output2=output;
% [~,index] = sortrows([output2.mean_FR].'); output2 = output2(index); clear index
% output2(1) = [];
% output2(78:79) = [];
clearvars -except output output2

for i = 1:length(output2)
    groups(i) = output2(i).age_DIV;
end
groups1 = unique(groups);

F1 = figure;
F1.OuterPosition = [750   533   456   467];

ind = find(groups ==groups1(1));

for i = 1:length(ind)
    data1(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(2));
for i = 1:length(ind)
    data2(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(3));
for i = 1:length(ind)
    data3(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(4));
for i = 1:length(ind)
    data4(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(5));
for i = 1:length(ind)
    data5(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(6));
for i = 1:length(ind)
    data6(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(7));
for i = 1:length(ind)
    data7(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(8));
for i = 1:length(ind)
    data8(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(9));
for i = 1:length(ind)
    data9(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(10));
for i = 1:length(ind)
    data10(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(11));
for i = 1:length(ind)
    data11(i) = output2(ind(i)).num_RC_nodes_mean;
end

clear ind; ind = find(groups ==groups1(12));
for i = 1:length(ind)
    data12(i) = output2(ind(i)).num_RC_nodes_mean;
end

notBoxPlot([data1],groups1(1),'style','sdline')
hold on
notBoxPlot([data2],groups1(2),'style','sdline')
notBoxPlot([data3],groups1(3),'style','sdline')
notBoxPlot([data4],groups1(4),'style','sdline')
notBoxPlot([data5],groups1(5),'style','sdline')
notBoxPlot([data6],groups1(6),'style','sdline')
notBoxPlot([data7],groups1(7),'style','sdline')
notBoxPlot([data8],groups1(8),'style','sdline')
notBoxPlot([data9],groups1(9),'style','sdline')
notBoxPlot([data10],groups1(10),'style','sdline')
notBoxPlot([data11],groups1(11),'style','sdline')
notBoxPlot([data12],groups1(12),'style','sdline')
% ylim([10 60])
aesthetics

%% firing rate ~ group * age

clear groups
for i = 1:length(output)
    
    groups{i} = strcat(num2str(output(i).age_DIV),output(i).grp);
    
end

group_names = unique(groups);

clear data
for j = 1:length(group_names)
    
    clear ind
    ind = strcmp(groups,group_names{j});
    ind = find(ind == 1);
    
    for i = 1:length(ind)
        data_temp(i,1) = output(ind(i)).mean_FR;
    end
    
    data{j} = data_temp;
    clear data_temp

end

% figure
% cols = [linspace(0,1,length(group_names))',linspace(1,0,length(group_names))',linspace(0,1,length(group_names))'];
% perm = [6     3     7     1     5    10     9    12     8     2    11     4];
% cols(:,1:3) = cols(perm,1:3);
% 
% rand
% for i = 1:length(group_names)
%     
%     grp = group_names{i};
%     b = notBoxPlot(data{i},str2num(grp(1:3)),'style','sdline');
%     b.data.MarkerFaceColor = cols(i,1:3);
%     hold on
%     clear grp
%     
% end
% aesthetics
% xlim([140 250])


figure
grp = group_names{1};
b = notBoxPlot(data{1},str2num(grp(1:3)),'style','sdline');
b.mu.Color = [1 1 1];
b.semPtch.FaceColor = [1 1 1]
b.sd.Color = [1 1 1];
b.semPtch.EdgeColor = [1 1 1];

hold on
b = notBoxPlot(data{2},str2num(grp(1:3))-0.2,'style','sdline');
b.data.Marker = 'd';
b.mu.Color = [1 1 1];
b.semPtch.FaceColor = [1 1 1]
b.sd.Color = [1 1 1];
b.semPtch.EdgeColor = [1 1 1];

hold on
b = notBoxPlot(data{3},str2num(grp(1:3))+0.2,'style','sdline');
b.data.Marker = 's';
b.mu.Color = [1 1 1];
b.semPtch.FaceColor = [1 1 1]
b.sd.Color = [1 1 1];
b.semPtch.EdgeColor = [1 1 1];



%% firing rate * SD * cell line
clear all
sd1 = load('manuel_3_FTDspikes3sd.mat');
sd2 = load('manuel_3_FTDspikes3.5sd.mat');
sd3 = load('manuel_3_FTDspikes4sd.mat');

clear sd1Gr
for i = 1:length(sd1.output)
    sd1FR(i) = sd1.output(i).mean_FR;
    sd1Gr(i) = sd1.output(i).grp;
    grouplist{i} = sd1.output(i).grp;
end

groups = unique(grouplist);
clear sd1FrGr
for i = 1:length(groups)
    sd1FrGr{:,i} = sd1FR(strfind(sd1Gr,groups{i}));    
end


clear sd2Gr grouplist
for i = 1:length(sd2.output)
    sd2FR(i) = sd2.output(i).mean_FR;
    sd2Gr(i) = sd2.output(i).grp;
    grouplist{i} = sd2.output(i).grp;
end

clear groups
groups = unique(grouplist);
clear sd2FrGr
for i = 1:length(groups)
    sd2FrGr{:,i} = sd2FR(strfind(sd2Gr,groups{i}));    
end


clear sd3Gr grouplist
for i = 1:length(sd3.output)
    sd3FR(i) = sd3.output(i).mean_FR;
    sd3Gr(i) = sd3.output(i).grp;
    grouplist{i} = sd3.output(i).grp;
end

clear groups
groups = unique(grouplist);
clear sd3FrGr
for i = 1:length(groups)
    sd3FrGr{:,i} = sd3FR(strfind(sd3Gr,groups{i}));    
end

f1 = figure;
hold on
for i = 1:length(groups)
scatter(i*ones(size(sd1FrGr{i})),sd1FrGr{i},10,'r')
end

for i = 1:length(groups)
scatter(i*ones(size(sd2FrGr{i})),sd2FrGr{i},10,'b')
end

for i = 1:length(groups)
scatter(i*ones(size(sd3FrGr{i})),sd3FrGr{i},10,'g')
end
aesthetics
xticks(1:4)
xticks(1:5)
xlabel('group')
ylabel('firing rate')


%% firing rate * SD * cell lines as separate lines
clear all; close all

sd1 = load('manuel_3_FTDspikes3sd.mat');
sd2 = load('manuel_3_FTDspikes3.5sd.mat');
sd3 = load('manuel_3_FTDspikes4sd.mat');

clear sd1Gr
for i = 1:length(sd1.output)
    sd1FR(i) = sd1.output(i).mean_FR;
    sd1Gr(i) = sd1.output(i).grp;
    grouplist{i} = sd1.output(i).grp;
end

groups = unique(grouplist);
clear sd1FrGr
for i = 1:length(groups)
    sd1FrGr{:,i} = sd1FR(strfind(sd1Gr,groups{i}));    
end


clear sd2Gr grouplist
for i = 1:length(sd2.output)
    sd2FR(i) = sd2.output(i).mean_FR;
    sd2Gr(i) = sd2.output(i).grp;
    grouplist{i} = sd2.output(i).grp;
end

clear groups
groups = unique(grouplist);
clear sd2FrGr
for i = 1:length(groups)
    sd2FrGr{:,i} = sd2FR(strfind(sd2Gr,groups{i}));    
end


clear sd3Gr grouplist
for i = 1:length(sd3.output)
    sd3FR(i) = sd3.output(i).mean_FR;
    sd3Gr(i) = sd3.output(i).grp;
    grouplist{i} = sd3.output(i).grp;
end

clear groups
groups = unique(grouplist);
clear sd3FrGr
for i = 1:length(groups)
    sd3FrGr{:,i} = sd3FR(strfind(sd3Gr,groups{i}));    
end

for i = 1:length(groups)
    scores1{i,1} = sd1FR(find(sd1Gr == groups{i}));
end

for i = 1:length(groups)
    scores2{i,1} = sd2FR(find(sd2Gr == groups{i}));
end

for i = 1:length(groups)
    scores3{i,1} = sd3FR(find(sd3Gr == groups{i}));
end

A = [scores1{1,1}; scores2{1,1}; scores3{1,1}];
B = [scores1{2,1}; scores2{2,1}; scores3{2,1}];
C = [scores1{3,1}; scores2{3,1}; scores3{3,1}];
D = [scores1{4,1}; scores2{4,1}; scores3{4,1}];
E = [scores1{5,1}; scores2{5,1}; scores3{5,1}];

f1 = figure;
colours = {[0.6350 0.0780 0.1840];[0.7 0.7 0.7];[0.2 0.2 0.2];	[1 0.5 0.5]};
hold on
for i  = 1:size(A,1)
    scatter((i*ones(1,length(A)))+0.15,A(i,:),50,colours{1,1},'o','jitter', 'off', 'jitterAmount', 0.3)
end

for i  = 1:size(B,1)
    scatter((i*ones(1,length(B)))+0.05,B(i,:),50,colours{2,1},'o','jitter', 'off', 'jitterAmount', 0.3)
end

for i  = 1:size(C,1)
    scatter((i*ones(1,length(C)))-0.05,C(i,:),50,colours{3,1},'o','jitter', 'off', 'jitterAmount', 0.3)
end

for i  = 1:size(D,1)
    scatter((i*ones(1,length(D)))-0.15,D(i,:),50,colours{4,1},'o','jitter', 'off', 'jitterAmount', 0.3)
end

% for i  = 1:size(E,1)
%     scatter(i*ones(1,length(E)),E(i,:),50,'y','o','jitter', 'on', 'jitterAmount', 0.3)
% end

pa = plot(-0.15  + [1 2 3], mean(A,2),'-s','MarkerFaceColor',colours{1,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{1,1});
pb = plot(-0.05 + [1 2 3], mean(B,2),'-s','MarkerFaceColor',colours{2,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{2,1});
pc = plot( 0.05 + [1 2 3], mean(C,2),'-s','MarkerFaceColor',colours{3,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{3,1});
pd = plot( 0.15  + [1 2 3], mean(D,2),'-s','MarkerFaceColor',colours{4,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{4,1});
% pe = plot(1.15  * [1 2 3], mean(E,2),'--sy','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10);
pa1 = plot(-0.15  + [1 2 3], mean(A,2),'s','MarkerFaceColor',colours{1,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{1,1});
pb1 = plot(-0.05 + [1 2 3], mean(B,2),'s','MarkerFaceColor',colours{2,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{2,1});
pc1 = plot( 0.05 + [1 2 3], mean(C,2),'s','MarkerFaceColor',colours{3,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{3,1});
pd1 = plot( 0.15  + [1 2 3], mean(D,2),'s','MarkerFaceColor',colours{4,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{4,1});
% legend([pa pb pc pd pe], ['CS30   '; 'EpiC   '; 'WTsli32'; 'CS29   '; 'CS30old'])
legend([pa pb pc pd], ['C9-L1'; 'H-L2 '; 'H-L1 '; 'C9-L2'],'box','off')


% std(A(1,:))) / (sqrt(length(A(1,:)))

aesthetics
xlim([0.5 3.5])
xticks(1:3)
xticklabels({'3.0';'3.5';'4.0'})
xlabel({'spike detection threshold','(std dev)'})
ylabel('firing rate (Hz)')
a = gca;
a.FontName = 'Arial';
a.FontSize = 18;

% saveas(f1,strcat('FR_bySD_byCellLine','.png'));
% saveas(f1,strcat('FR_bySD_byCellLine','.epsc'));

FR_30    = [A(1,:) B(1,:) C(1,:) D(1,:)]';
FR_35    = [A(2,:) B(2,:) C(2,:) D(2,:)]';
FR_40    = [A(3,:) B(3,:) C(3,:) D(3,:)]';
grp = erase(sd1Gr,'E');

% one-way anova for each sync window
[p,tbl,stats] = anova1(FR_30,grp'); close(gcf); close(gcf)
Fstat = tbl{2,5};
% add to plot anova p value
text(0.75,25,strcat({'p = '},num2str(round(p,3))),'FontSize',12)
% post hoc tukey
% multcompare(stats)% last column is p value, first 2 columns indicate
% which groups are being compared
[p,tbl,stats] = anova1(FR_35,grp'); close(gcf); close(gcf)
Fstat = tbl{2,5};
% add to plot anova p value
text(1.75,14,strcat({'p = '},num2str(round(p,3))),'FontSize',12)

[p,tbl,stats] = anova1(FR_40,grp'); close(gcf); close(gcf)
Fstat = tbl{2,5};
% add to plot anova p value
text(2.75,4,strcat({'p = '},num2str(round(p,3))),'FontSize',12)

saveas(f1,strcat('FR_bySD_byCellLine','.png'));
saveas(f1,strcat('FR_bySD_byCellLine','.epsc'));

%% STTC * STTC_synthetic * cell lines as separate lines
clear all
files = dir('*FTD*mSpikes_3_CTRL_adjM_0.175.mat');
active_only = 1; %set to one to include on channels with connections

for i = 1:length(files)
    load(files(i).name)
    filename = files(i).name;
    adjMs = adjM2s;
%     % run this loop if self connections still exist in the matrix
%     for j = 1:size (adjM2s,3) 
%         adjMs(:,:,j) = adjM2s(:,:,j) - eye(size(adjM2s(:,:,j)));
%     end   
    if active_only == 1
        active_chan_index = find(sum(adjMs(:,:,1)) ~= 0);
        adjMs = adjMs(active_chan_index,active_chan_index,:);
    end
%     % loop to check it's the same for all ctrl matrices
%     % (it should be)
%     for k=1:size(w,3) 
%         active_chan_indices(k,:) = find(w(:,:,k) ~= 0);
%     end 
%     w = mean(adjMs); %%this way does not exclude self-connection value from contributing to mean
    w = sum(adjMs) ./ (size(adjMs,1)-1) ; %this way excludes self-connection value from contributing to mean
    w(find(isnan(w) == 1)) = 0; % remove NaNs
    w_synt(i) = mean( reshape(w, [size(w,1), size(w,2) * size(w,3)]) );
     
    strRemove = '_CTRL';
    real_file = erase(filename,strRemove);
    
    %find group info
    if ~isempty(strfind(real_file(1:end-4),'Grp'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(real_file(1:end-4),'Group'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Group')+5);
    else
        grp (i) = ['E'];
    end
    
    load(real_file)
    adjM = adjM - eye(size(adjM));
    adjM(find(isnan(adjM) == 1)) = 0; % remove NaNs
    if active_only == 1
        active_chan_index = find(sum(adjM) ~= 0);
        adjM = adjM(active_chan_index,active_chan_index);
    end
    w_real(i) = mean(  sum(adjM) ./ (length(adjM)-1)  );          
end

clearvars -except w_real w_synt grp

pcnt_sig = 100 * (length(find(w_real > prctile(w_synt,95))) / length(w_real));
% i.e. percent that occurred only 5% or less of the time in synthetic
% trains
grouplist = grp;

groups = unique(grouplist);

for i = 1:length(groups)
    scores1{i,1} = w_real(find(grp == groups(i)));
end

for i = 1:length(groups)
    scores2{i,1} = w_synt(find(grp == groups(i)));
end


A = [scores1{1,1}; scores2{1,1}];
B = [scores1{2,1}; scores2{2,1}];
C = [scores1{3,1}; scores2{3,1}];
D = [scores1{4,1}; scores2{4,1}];
E = [scores1{5,1}; scores2{5,1}];

f1 = figure;
colours = {[0.6350 0.0780 0.1840];[0.7 0.7 0.7];[0.2 0.2 0.2];	[1 0.5 0.5]};
hold on
for i  = 1:size(A,1)
    scatter((i*ones(1,length(A)))-0.15,A(i,:),50,colours{1,1},'o','jitter', 'off', 'jitterAmount', 0.1,...
        'LineWidth',1.5,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4)
end

for i  = 1:size(B,1)
    scatter((i*ones(1,length(B)))-0.05,B(i,:),50,colours{2,1},'o','jitter', 'off', 'jitterAmount', 0.1,...
        'LineWidth',1.5,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4)
end

for i  = 1:size(C,1)
    scatter((i*ones(1,length(C)))+0.05,C(i,:),50,colours{3,1},'o','jitter', 'off', 'jitterAmount', 0.1,...
        'LineWidth',1.5,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4)
end

for i  = 1:size(D,1)
    scatter((i*ones(1,length(D)))+0.15,D(i,:),50,colours{4,1},'o','jitter', 'off', 'jitterAmount', 0.1,...
        'LineWidth',1.5,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4)
end

% for i  = 1:size(E,1)
%     scatter(i*ones(1,length(E)),E(i,:),50,'y','o','jitter', 'on', 'jitterAmount', 0.3)
% end

% error bars
sema = [std(A(1,:)) / sqrt(length(A)) , std(A(2,:)) / sqrt(length(A))];
ea = errorbar(-0.15  + [1 2],mean(A,2), sema,'LineStyle','none','Color','k','LineWidth',2);
semb = [std(B(1,:)) / sqrt(length(B)) , std(B(2,:)) / sqrt(length(B))];
eb = errorbar(-0.05 + [1 2],mean(B,2), sema,'LineStyle','none','Color','k','LineWidth',2);
semc = [std(C(1,:)) / sqrt(length(C)) , std(C(2,:)) / sqrt(length(C))];
ec = errorbar( 0.05 + [1 2],mean(C,2), sema,'LineStyle','none','Color','k','LineWidth',2);
semd = [std(D(1,:)) / sqrt(length(D)) , std(D(2,:)) / sqrt(length(D))];
ed = errorbar( 0.15  + [1 2],mean(D,2), sema,'LineStyle','none','Color','k','LineWidth',2);
% seme = [std(E(1,:)) / sqrt(length(E)) , std(E(2,:)) / sqrt(length(E))];
% ee = errorbar(1.15  * [1 2 3],mean(E,2), sema,'LineStyle','none','Color','k','LineWidth',2);

% pa = plot(-0.15  + [1 2], mean(A,2),'--sr','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2);
% pb = plot(-0.05 + [1 2], mean(B,2),'--sc','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2);
% pc = plot( 0.05 + [1 2], mean(C,2),'--sb','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2);
% pd = plot( 0.15  + [1 2], mean(D,2),'--sm','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2);
% pe = plot(1.15  * [1 2 3], mean(E,2),'--sy','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2);

% legend([pa pb pc pd pe], ['CS30   '; 'EpiC   '; 'WTsli32'; 'CS29   '; 'CS30old'])
% legend([pa pb pc pd], ['CS30   '; 'EpiC   '; 'WTsli42'; 'CS29   '],'box','off')
pa = plot(-0.15  + [1 2], mean(A,2),'-s','MarkerFaceColor',colours{1,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{1,1});
pb = plot(-0.05 + [1 2], mean(B,2),'-s','MarkerFaceColor',colours{2,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{2,1});
pc = plot( 0.05 + [1 2], mean(C,2),'-s','MarkerFaceColor',colours{3,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{3,1});
pd = plot( 0.15  + [1 2], mean(D,2),'-s','MarkerFaceColor',colours{4,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{4,1});
% pe = plot(1.15  * [1 2 3], mean(E,2),'--sy','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10);
pa1 = plot(-0.15  + [1 2], mean(A,2),'s','MarkerFaceColor',colours{1,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{1,1});
pb1 = plot(-0.05 + [1 2], mean(B,2),'s','MarkerFaceColor',colours{2,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{2,1});
pc1 = plot( 0.05 + [1 2], mean(C,2),'s','MarkerFaceColor',colours{3,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{3,1});
pd1 = plot( 0.15  + [1 2], mean(D,2),'s','MarkerFaceColor',colours{4,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{4,1});

% std(A(1,:))) / (sqrt(length(A(1,:)))

aesthetics
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'real';'synth.'})
% xlabel('multiplier')
ylabel('global functional connectivity')
a = gca;
a.FontName = 'Arial';
a.FontSize = 18;

% add p values
bwGroupsANOVA_real_p = 0.226; % from R analysis
mainEffectSource     = 0.001; % from R analysis
x = [0.8:0.1:1.2]';
y = 0.38 * ones(length(x),1);
plot(x,y,'-k','LineWidth',1); text(x(1),0.39,strcat({'p = '},num2str(bwGroupsANOVA_real_p)),'FontSize',12);
x = [1 2]';
y = 0.40 * ones(length(x),1);
plot(x,y,'-k','LineWidth',1); text(x(1),0.41,strcat({'p < '},num2str(mainEffectSource)),'FontSize',12);

% legend([pa pb pc pd pe], ['CS30   '; 'EpiC   '; 'WTsli32'; 'CS29   '; 'CS30old'])
legend([pa pb pc pd], ['C9-L1'; 'H-L2 '; 'H-L1 '; 'C9-L2'],'box','off')

saveas(f1,strcat('STTC_byType_byCellLine','.png'));
saveas(f1,strcat('STTC_byType_byCellLine','.epsc'));


%% STTC * cell lines; as bars instead of squares for means
% (no synthetic data for now, could later add in dashed line the width of
% each bar as the mean of the synthetics to show chance level
clear all
files = dir('*FTD*mSpikes_3_CTRL_adjM_0.175.mat');
active_only = 1; %set to one to include on channels with connections

for i = 1:length(files)
    load(files(i).name)
    filename = files(i).name;
    adjMs = adjM2s;
%     % run this loop if self connections still exist in the matrix
%     for j = 1:size (adjM2s,3) 
%         adjMs(:,:,j) = adjM2s(:,:,j) - eye(size(adjM2s(:,:,j)));
%     end   
    if active_only == 1
        active_chan_index = find(sum(adjMs(:,:,1)) ~= 0);
        adjMs = adjMs(active_chan_index,active_chan_index,:);
    end
%     % loop to check it's the same for all ctrl matrices
%     % (it should be)
%     for k=1:size(w,3) 
%         active_chan_indices(k,:) = find(w(:,:,k) ~= 0);
%     end 
%     w = mean(adjMs); %%this way does not exclude self-connection value from contributing to mean
    w = sum(adjMs) ./ (size(adjMs,1)-1) ; %this way excludes self-connection value from contributing to mean
    w(find(isnan(w) == 1)) = 0; % remove NaNs
    w_synt(i) = mean( reshape(w, [size(w,1), size(w,2) * size(w,3)]) );
     
    strRemove = '_CTRL';
    real_file = erase(filename,strRemove);
    
    %find group info
    if ~isempty(strfind(real_file(1:end-4),'Grp'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(real_file(1:end-4),'Group'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Group')+5);
    else
        grp (i) = ['E'];
    end
    
    load(real_file)
    adjM = adjM - eye(size(adjM));
    adjM(find(isnan(adjM) == 1)) = 0; % remove NaNs
    if active_only == 1
        active_chan_index = find(sum(adjM) ~= 0);
        adjM = adjM(active_chan_index,active_chan_index);
    end
    w_real(i) = mean(  sum(adjM) ./ (length(adjM)-1)  );          
end

clearvars -except w_real w_synt grp

pcnt_sig = 100 * (length(find(w_real > prctile(w_synt,95))) / length(w_real));
% i.e. percent that occurred only 5% or less of the time in synthetic
% networks
grouplist = grp;

groups = unique(grouplist);

for i = 1:length(groups)
    scores1{i,1} = w_real(find(grp == groups(i)));
end

for i = 1:length(groups)
    scores2{i,1} = w_synt(find(grp == groups(i)));
end


A = [scores1{1,1}; scores2{1,1}];
B = [scores1{2,1}; scores2{2,1}];
C = [scores1{3,1}; scores2{3,1}];
D = [scores1{4,1}; scores2{4,1}];
E = [scores1{5,1}; scores2{5,1}];

f1 = figure;
hold on

bar_scores = [mean(C(1,:)) mean(B(1,:)) mean(A(1,:)) mean(D(1,:))];
b = bar(1:length(bar_scores), bar_scores,0.5);
% re color bars
ECvec = {0.5*[1 1 1],0.5*[1 1 1],[1 0 0],[1 0 0]};
hold on
jitval = 0.625;
for  i=1:length(bar_scores)
    bar(i, bar_scores(i),jitval,'FaceColor','w', 'EdgeColor',ECvec{i},'LineWidth',1.5);
end
xlim([0.5 4.5])
xticks(1:length(bar_scores))
ylim(0.1 *   [1   ceil(max(bar_scores)*10)])
yticks(0.1 * [1 : ceil(max(bar_scores)*10)])
aesthetics
% % re color axis
plot(xlim,[min(ylim) min(ylim)],'-k','LineWidth',2)
% bar(1:length(bar_scores), bar_scores,0.5,'Visible','off');

bx = gca;
bx.TickDir = 'out';
bx.TickLength = 0.02 * [1 1];
bx.LineWidth = 2;
bx.FontName = 'Arial';
bx.FontSize = 18;
f1.Position = [1030 270 424 420];
ylabel('global connectivity')
xticklabels({'WT  ','EpiC','CS30','CS29'})
xtickangle(45)

% make sure that: 
% C is WT42, B is EpiC, A is CS30, D is CS29
scat_scores = {C(1,:); B(1,:); A(1,:); D(1,:)};
FCvec = {0.5*[1 1 1],[1 1 1],[1 0 0],[1 1 1]};
for i = 1:length(scat_scores)
    scatter(i*ones(1,length(scat_scores{i})),scat_scores{i,:},50,ECvec{i},'o',...
        'jitter', 'on', 'jitterAmount', jitval/2.5,...
        'LineWidth',1.5,'MarkerFaceColor',FCvec{i},'MarkerEdgeColor',ECvec{i})
end

% error bars
for i = 1:length(scat_scores)
    SEMs(i) = [std(scat_scores{i})   /   (sqrt(length(scat_scores{i})))]; 
end

SEMcolVec = ECvec;
for i = 1:length(scat_scores)
%     errorbar(i,mean(scat_scores{i}),SEMs(i),'LineStyle','none','Color',ECvec{i},'LineWidth',1.5,'CapSize',18)
    errorbar(i,mean(scat_scores{i}),SEMs(i),'LineStyle','none','Color','k','LineWidth',1.5,'CapSize',18)
end

saveas(f1,strcat('STTC_Bar_byCellLine','.png'));
saveas(f1,strcat('STTC_Bar_byCellLine','.epsc'));

%% supplementary table info
clear all
files = dir('*FTD*mSpikes_3_CTRL_adjM_0.175.mat');
active_only = 1; %set to one to include on channels with connections


for i = 1:length(files)
    load(files(i).name)
    filename = files(i).name;
    adjMs = adjM2s;
    %     % run this loop if self connections still exist in the matrix
    %     for j = 1:size (adjM2s,3)
    %         adjMs(:,:,j) = adjM2s(:,:,j) - eye(size(adjM2s(:,:,j)));
    %     end
    if active_only == 1
        active_chan_index = find(sum(adjMs(:,:,1)) ~= 0);
        adjMs = adjMs(active_chan_index,active_chan_index,:);
    end
    %     % loop to check it's the same for all ctrl matrices
    %     % (it should be)
    %     for k=1:size(w,3)
    %         active_chan_indices(k,:) = find(w(:,:,k) ~= 0);
    %     end
    %     w = mean(adjMs); %%this way does not exclude self-connection value from contributing to mean
    w = sum(adjMs) ./ (size(adjMs,1)-1) ; %this way excludes self-connection value from contributing to mean
    w(find(isnan(w) == 1)) = 0; % remove NaNs
    
    % get all 100 synthetic networks' global connectivity for each rec. to
    % calculate p values for real networks (i.e. permutation test)
    % w is the mean weight for each channel across 100 iterations
    gc = mean(w);
    glob_con(i,:) = reshape(gc, [size(gc,1), size(gc,2) * size(gc,3)]);
    
    w_synt(i) = mean( reshape(w, [size(w,1), size(w,2) * size(w,3)]) );
    
    strRemove = '_CTRL';
    real_file = erase(filename,strRemove);
    
    %find group info
    if ~isempty(strfind(real_file(1:end-4),'Grp'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(real_file(1:end-4),'Group'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Group')+5);
    else
        grp (i) = ['E'];
    end
    
    load(real_file)
    adjM = adjM - eye(size(adjM));
    adjM(find(isnan(adjM) == 1)) = 0; % remove NaNs
    if active_only == 1
        active_chan_index = find(sum(adjM) ~= 0);
        adjM = adjM(active_chan_index,active_chan_index);
    end
    w_real(i) = mean(  sum(adjM) ./ (length(adjM)-1)  );
    
    %get p value for Glob connectivity for each recording
    sorted_gc = sort(glob_con(i,:));
    p_val(i)  = (length(find(sorted_gc > w_real(i))) )  /  length(sorted_gc);
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get n nodes with degree > k
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    count2 = 1;     
    
    for cutoff = [60 75 90]
        
        threshold = prctile(adjM(:), cutoff);
        
        edges=adjM;
        edges(find(edges >= threshold))= 1;
        edges(find(edges < threshold)) = 0;
        
        DegreeVec(:,count2)       = sum(edges);
        
        edges_synt = adjMs;
        edges_synt(find(adjMs >= threshold)) = 1;
        edges_synt(find(adjMs <  threshold)) = 0;
        
        DegreeVec_synt_all(:,:,count2)  = sum(edges_synt);        
        
        count2 = count2 + 1;
        
    end
    
    dg{i,1}         = round(mean(DegreeVec,2));
    dg_synt{i,1}    = round(mean(DegreeVec_synt_all,3)); % gets the mean across all threshold for each channel for each iteration of permutation
    clear DegreeVec edges threshold DegreeVec_synt_all
    
end

clearvars -except w_real w_synt grp  p_val dg files dg_synt 




    %%%%%%%%%%%%%%%%%%%%
    % num nodes with k > 30
    %%%%%%%%%%%%%%%%%%%%
    for i = 1 : length(files)
        nk30_real(i) = length(find(dg{i,1} >= 30));
        nk30_synt(i) = length(find(dg_synt{i,1} >= 30))  /  length(dg_synt{i,1});
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  get m and sem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_vec  = grp(find(p_val <= 0.05));
sig_vec2 = grp(find(p_val >  0.05));
grouplist = grp;
groups = unique(grouplist);
for i = 1 : length(groups)
    n(i) = length(find(grp == groups(i))); 
    percent_sig(i)      =  round( 100*(length( find (sig_vec == groups(i) ) ) )  /  (n(i)) , 1); 
    mean_pval(i)        = mean(p_val(find(groups(i) == grp)));
    sem_pval(i)         = (std(p_val(find(groups(i) == grp))))     / (sqrt(n(i)))   ;
    
    mean_GC(i)          = mean(w_real(find(groups(i) == grp)));
    sem_GC(i)           = (std(w_real(find(groups(i) == grp))))     / (sqrt(n(i)))   ;
    
    mean_GCsynt(i)      = mean(w_synt(find(groups(i) == grp)));
    sem_GCsynt(i)       = (std(w_synt(find(groups(i) == grp))))     / (sqrt(n(i)))   ;
    
    mean_nk30(i)        =  mean(nk30_real(find(groups(i) == grp)));
    sem_nk30(i)         = (std(nk30_real(find(groups(i) == grp))))     / (sqrt(n(i)))   ;
    
    mean_nk30synt(i)    =  mean(nk30_synt(find(groups(i) == grp)));
    sem_nk30synt(i)     = (std(nk30_synt(find(groups(i) == grp))))     / (sqrt(n(i)))   ;

end

% note: A = CS30  |  B = EpiC  |  C = WTsli42  |  D = CS29  |  E = CS30 older

% for i = 1 : length(grp)
%     
%     cellgrp{i,1} = grp(i);
%     
%     if strcmp(cellgrp{i,1} , 'A')
%         cellgrp{i,2} = 'CS30';
%     elseif strcmp(cellgrp{i,1} , 'B')
%         cellgrp{i,2} = 'EpiC';
%     elseif strcmp(cellgrp{i,1} , 'C')
%         cellgrp{i,2} = 'WTsli42';
%     elseif strcmp(cellgrp{i,1} , 'D')
%         cellgrp{i,2} = 'CS29';
%     elseif strcmp(cellgrp{i,1} , 'E')
%         cellgrp{i,2} = 'CS30older';
%     end
%     
% end
% 
% tabgrp = cell2table(cellgrp);
clear grpstruct tabgrp2
for i = 1 : length(grp)
    
    grpstruct(i).file = files(i).name;
    grpstruct(i).line = grp(i);
    
    if strcmp(grpstruct(i).line , 'A')
        grpstruct(i).label = "CS30";
    elseif strcmp(grpstruct(i).line , 'B')
        grpstruct(i).label = "EpiC";
    elseif strcmp(grpstruct(i).line , 'C')
        grpstruct(i).label = "WTsli42";
    elseif strcmp(grpstruct(i).line , 'D')
        grpstruct(i).label = "CS29";
    elseif strcmp(grpstruct(i).line , 'E')
        grpstruct(i).label = "CS30older";
    end
    
end

tabgrp2 = struct2table(grpstruct);
% writetable(tabgrp2, strcat('FTD_groups_files','.xlsx'))


% for i = 1:length(files)
%    suppstats(i).line = grp(i);
%    suppstats(i).
% end

vals(:,1) = nk30_real';
vals(:,2) = nk30_synt';
vals(:,3) = w_real';
vals(:,4) = w_synt';
% writematrix(vals, strcat('suppstats','.xlsx'))




%%%%%%%%%%%%%%% plot fig for n k>30
for i = 1:length(groups)
    scores1{i,1} = nk30_real(find(grp == groups(i)));
end

for i = 1:length(groups)
    scores2{i,1} = nk30_synt(find(grp == groups(i)));
end


A = [scores1{1,1}; scores2{1,1}];
B = [scores1{2,1}; scores2{2,1}];
C = [scores1{3,1}; scores2{3,1}];
D = [scores1{4,1}; scores2{4,1}];
E = [scores1{5,1}; scores2{5,1}];

f1 = figure;
colours = {[0.6350 0.0780 0.1840];[0.7 0.7 0.7];[0.2 0.2 0.2];	[1 0.5 0.5]};
hold on

for i  = 1:size(A,1)
    scatter((i*ones(1,length(A)) )-0.15,A(i,:),50,colours{1,1},'o','jitter', 'off', 'jitterAmount', 0.1,...
        'LineWidth',1.5,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4)
end

for i  = 1:size(B,1)
    scatter((i*ones(1,length(B)) )-0.05,B(i,:),50,colours{2,1},'o','jitter', 'off', 'jitterAmount', 0.1,...
        'LineWidth',1.5,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4)
end

for i  = 1:size(C,1)
    scatter((i*ones(1,length(C)) )+0.05,C(i,:),50,colours{3,1},'o','jitter', 'off', 'jitterAmount', 0.1,...
        'LineWidth',1.5,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4)
end

for i  = 1:size(D,1)
    scatter((i*ones(1,length(D)) )+0.15,D(i,:),50,colours{4,1},'o','jitter', 'off', 'jitterAmount', 0.1,...
        'LineWidth',1.5,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4)
end

% for i  = 1:size(E,1)
%     scatter(i*ones(1,length(E)),E(i,:),50,'y','o','jitter', 'on', 'jitterAmount', 0.3)
% end

% error bars
sema = [std(A(1,:)) / sqrt(length(A)) , std(A(2,:)) / sqrt(length(A))];
ea = errorbar(-0.15  + [1 2],mean(A,2), sema,'LineStyle','none','Color','k','LineWidth',2);
semb = [std(B(1,:)) / sqrt(length(B)) , std(B(2,:)) / sqrt(length(B))];
eb = errorbar(-0.05 + [1 2],mean(B,2), sema,'LineStyle','none','Color','k','LineWidth',2);
semc = [std(C(1,:)) / sqrt(length(C)) , std(C(2,:)) / sqrt(length(C))];
ec = errorbar( 0.05 + [1 2],mean(C,2), sema,'LineStyle','none','Color','k','LineWidth',2);
semd = [std(D(1,:)) / sqrt(length(D)) , std(D(2,:)) / sqrt(length(D))];
ed = errorbar( 0.15  + [1 2],mean(D,2), sema,'LineStyle','none','Color','k','LineWidth',2);
% seme = [std(E(1,:)) / sqrt(length(E)) , std(E(2,:)) / sqrt(length(E))];
% ee = errorbar(1.15  * [1 2 3],mean(E,2), sema,'LineStyle','none','Color','k','LineWidth',2);

pa = plot(-0.15  + [1 2], mean(A,2),'-s','MarkerFaceColor',colours{1,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{1,1});
pb = plot(-0.05 + [1 2], mean(B,2),'-s','MarkerFaceColor',colours{2,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{2,1});
pc = plot( 0.05 + [1 2], mean(C,2),'-s','MarkerFaceColor',colours{3,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{3,1});
pd = plot( 0.15  + [1 2], mean(D,2),'-s','MarkerFaceColor',colours{4,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{4,1});
% pe = plot(1.15  * [1 2 3], mean(E,2),'--sy','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10);
pa1 = plot(-0.15  + [1 2], mean(A,2),'s','MarkerFaceColor',colours{1,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{1,1});
pb1 = plot(-0.05 + [1 2], mean(B,2),'s','MarkerFaceColor',colours{2,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{2,1});
pc1 = plot( 0.05 + [1 2], mean(C,2),'s','MarkerFaceColor',colours{3,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{3,1});
pd1 = plot( 0.15  + [1 2], mean(D,2),'s','MarkerFaceColor',colours{4,1},...
    'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'Color',colours{4,1});
% legend([pa pb pc pd pe], ['CS30   '; 'EpiC   '; 'WTsli32'; 'CS29   '; 'CS30old'])

aesthetics
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'real';'synth.'})
% xlabel('multiplier')
ylabel('# high-degree nodes')
a = gca;
a.FontName = 'Arial';
a.FontSize = 18;

% add p values
bwGroupsANOVA_real_p = 0.110; % from R analysis
% mainEffectSource     = 0.003; % from R analysis; doesn't treat source as rep. measures
mainEffectSource     = 0.008; % from SPSS analysis; treats source as rep. measures
x = [0.8:0.1:1.2]';
y = 9.3 * ones(length(x),1);
plot(x,y,'-k','LineWidth',1); text(x(1),9.65,strcat({'p = '},num2str(bwGroupsANOVA_real_p)),'FontSize',12);
x = [1 2]';
y = 10 * ones(length(x),1);
plot(x,y,'-k','LineWidth',1); text(x(1),10.35,strcat({'p = '},num2str(mainEffectSource)),'FontSize',12);

l=legend([pa pb pc pd], ['C9-L1'; 'H-L2 '; 'H-L1 '; 'C9-L2'],'box','off');
l.Location = 'north';
% std(A(1,:))) / (sqrt(length(A(1,:)))

saveas(f1,strcat('nk30_byType_byCellLine','.png'));
saveas(f1,strcat('nk30_byType_byCellLine','.epsc'));

%% sttc over cell line and sync window

clear all
files = dir('*FTD*mSpikes_3_adjM_0.125.mat'); files = files(~contains({files.name}, '190814'));
active_only = 0;

for i = 1:length(files)
    load(files(i).name)
    real_file = files(i).name;
    
    %find group info
    if ~isempty(strfind(real_file(1:end-4),'Grp'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(real_file(1:end-4),'Group'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Group')+5);
    else
        grp (i) = ['E'];
    end
    
    load(real_file)
    adjM = adjM - eye(size(adjM));
    adjM(find(isnan(adjM) == 1)) = 0; % remove NaNs
    if active_only == 1
        active_chan_index = find(sum(adjM) ~= 0);
        adjM = adjM(active_chan_index,active_chan_index);
    end
    w_125(i) = mean(  sum(adjM) ./ (length(adjM)-1)  );
end

clearvars -except w_125 active_only
files = dir('*FTD*mSpikes_3_adjM_0.15.mat'); files = files(~contains({files.name}, '190814'));

for i = 1:length(files)
    load(files(i).name)
    real_file = files(i).name;
    
    %find group info
    if ~isempty(strfind(real_file(1:end-4),'Grp'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(real_file(1:end-4),'Group'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Group')+5);
    else
        grp (i) = ['E'];
    end
    
    load(real_file)
    adjM = adjM - eye(size(adjM));
    adjM(find(isnan(adjM) == 1)) = 0; % remove NaNs
    if active_only == 1
        active_chan_index = find(sum(adjM) ~= 0);
        adjM = adjM(active_chan_index,active_chan_index);
    end
    w_150(i) = mean(  sum(adjM) ./ (length(adjM)-1)  );
end

clearvars -except w_125 w_150 active_only
files = dir('*FTD*mSpikes_3_adjM_0.175.mat'); files = files(~contains({files.name}, '190814'));

for i = 1:length(files)
    load(files(i).name)
    real_file = files(i).name;
    
    %find group info
    if ~isempty(strfind(real_file(1:end-4),'Grp'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(real_file(1:end-4),'Group'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Group')+5);
    else
        grp (i) = ['E'];
    end
    
    load(real_file)
    adjM = adjM - eye(size(adjM));
    adjM(find(isnan(adjM) == 1)) = 0; % remove NaNs
    if active_only == 1
        active_chan_index = find(sum(adjM) ~= 0);
        adjM = adjM(active_chan_index,active_chan_index);
    end
    w_175(i) = mean(  sum(adjM) ./ (length(adjM)-1)  );
end

clearvars -except w_125 w_150 w_175 active_only
files = dir('*FTD*mSpikes_3_adjM_0.2.mat'); files = files(~contains({files.name}, '190814'));

for i = 1:length(files)
    load(files(i).name)
    real_file = files(i).name;
    
    %find group info
    if ~isempty(strfind(real_file(1:end-4),'Grp'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(real_file(1:end-4),'Group'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Group')+5);
    else
        grp (i) = ['E'];
    end
    
    load(real_file)
    adjM = adjM - eye(size(adjM));
    adjM(find(isnan(adjM) == 1)) = 0; % remove NaNs
    if active_only == 1
        active_chan_index = find(sum(adjM) ~= 0);
        adjM = adjM(active_chan_index,active_chan_index);
    end
    w_200(i) = mean(  sum(adjM) ./ (length(adjM)-1)  );
end

clearvars -except w_125 w_150 w_175 w_200 files grp


%%%%%%%%%%%%%%%%%% get groups
grouplist = grp;
groups = unique(grouplist);

w = [w_125' w_150' w_175' w_200'];



for i = 1 : length(groups)
        
        n(i) = length(find(grp == groups(i)));
        
        group_means(i,:) = mean(w(find(grp == groups(i)),:));
        group_SEMs(i,:)  = (std(w(find(grp == groups(i)),:)))     / (sqrt(n(i)))   ;
               
end





f1 = figure;
hold on
timewins = 1:size(w,2);
% error bars
% for timwin = 1 : size(w,2)
stagger = [-0.15 -0.05 0.05 0.15];
colours = {[0.6350 0.0780 0.1840];[0.7 0.7 0.7];[0.2 0.2 0.2];	[1 0.5 0.5]};
    
for group = 1 : length(groups)
    
    ea = errorbar(stagger(group)  + timewins,group_means(group,:), ...
        group_SEMs(group,:),...
        'LineStyle','none','Color','k','LineWidth',2);
    
    pa(i) = plot(stagger(group)  + timewins,group_means(group,:),...
        '-s','Color',colours{group,1},'MarkerFaceColor',colours{group,1},'MarkerEdgeColor','k',...
        'MarkerSize',10,'LineWidth',2);
    
end
    
for group = 1 : length(groups)
    
    ea = errorbar(stagger(group)  + timewins,group_means(group,:), ...
        group_SEMs(group,:),...
        'LineStyle','none','Color','k','LineWidth',2);
    
    pa(i) = plot(stagger(group)  + timewins,group_means(group,:),...
        's','Color',colours{group,1},'MarkerFaceColor',colours{group,1},'MarkerEdgeColor','k',...
        'MarkerSize',10,'LineWidth',2);
    
end
    
% legend([pa pb pc pd pe], ['CS30   '; 'EpiC   '; 'WTsli32'; 'CS29   '; 'CS30old'])
L1 = plot(10  + timewins,group_means(1,:),...
        '-s','Color',colours{1,1},'MarkerFaceColor',colours{1,1},'MarkerEdgeColor','k',...
        'MarkerSize',10,'LineWidth',2,'Visible','on');
L2 = plot(10  + timewins,group_means(2,:),...
        '-s','Color',colours{2,1},'MarkerFaceColor',colours{2,1},'MarkerEdgeColor','k',...
        'MarkerSize',10,'LineWidth',2,'Visible','on');
L3 = plot(10  + timewins,group_means(3,:),...
        '-s','Color',colours{3,1},'MarkerFaceColor',colours{3,1},'MarkerEdgeColor','k',...
        'MarkerSize',10,'LineWidth',2,'Visible','on');
L4 = plot(10  + timewins,group_means(4,:),...
        '-s','Color',colours{4,1},'MarkerFaceColor',colours{4,1},'MarkerEdgeColor','k',...
        'MarkerSize',10,'LineWidth',2,'Visible','on');

leg1 = legend([L1 L2 L3 L4],['C9-L1'; 'H-L2 '; 'H-L1 '; 'C9-L2'],'Location','northwest',...
    'box','off');
% set(leg1, 'SelectionHighlight', 'off');
% https://uk.mathworks.com/matlabcentral/answers/233998-set-legend-to-non-transparent
% std(A(1,:))) / (sqrt(length(A(1,:)))

aesthetics
xlim([0.5 4.5])
xticks(1:4)
xticklabels({'125';'150';'175';'200'})
xlabel('synchronicity window (ms)')
ylabel('global functional connectivity')
a = gca;
a.FontName = 'Arial';
a.FontSize = 16;

% saveas(f1,strcat('global_fc_byWindow_byCellLine','.png'));
% saveas(f1,strcat('global_fc_byWindow_byCellLine','.epsc'));

% one-way anova for each sync window
[p,tbl,stats] = anova1(w_125,grp'); close(gcf); close(gcf)
Fstat = tbl{2,5};
% add to plot anova p value
text(0.75,max(w_125)+0.15,strcat({'p = '},num2str(round(p,3))),'FontSize',12)
% post hoc tukey
% multcompare(stats)% last column is p value, first 2 columns indicate
% which groups are being compared
[p,tbl,stats] = anova1(w_150,grp'); close(gcf); close(gcf)
Fstat = tbl{2,5};
% add to plot anova p value
text(1.75,max(w_150)+0.15,strcat({'p = '},num2str(round(p,3))),'FontSize',12)

[p,tbl,stats] = anova1(w_175,grp'); close(gcf); close(gcf)
Fstat = tbl{2,5};
% add to plot anova p value
text(2.75,max(w_175)+0.15,strcat({'p = '},num2str(round(p,3))),'FontSize',12)

[p,tbl,stats] = anova1(w_200,grp'); close(gcf); close(gcf)
Fstat = tbl{2,5};
% add to plot anova p value
text(3.75,0.65,strcat({'p = '},num2str(round(p,3))),'FontSize',12)
% ylim([0 0.8])

saveas(f1,strcat('global_fc_byWindow_byCellLine','.png'));
saveas(f1,strcat('global_fc_byWindow_byCellLine','.epsc'));

%% supplementary table FIRING RATE
clear all
files = dir('*FTD*mSpikes_3_CTRL_spikeMat.mat');
active_only = 1; %set to one to include on channels with connections
rec_length_s = 360;

progressbar('files'); p = gcf;

for i = 1:length(files)
    load(files(i).name)
    filename = files(i).name;
    
    FR_synt(i) = mean( full(sum(rand_spikeMat_sparse_all{1,1})) ./ rec_length_s );
    
    clear rand_spikeMat_sparse_all
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % real file
    
    strRemove = '_CTRL_spikeMat';
    real_file = erase(filename,strRemove);
    
    %find group info
    if ~isempty(strfind(real_file(1:end-4),'Grp'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(real_file(1:end-4),'Group'))
        grp (i) = real_file(strfind(real_file(1:end-4),'Group')+5);
    else
        grp (i) = ['E'];
    end
    
    load(real_file)
    
    FR_real(i) = mean( full(sum(mSpikes)) ./ rec_length_s );
    
    clear filename real_file mSpikes
    
    progressbar(i/length(files))
    
end

clearvars -except FR_real FR_synt grp files  


grouplist = grp;
groups = unique(grouplist);
for i = 1 : length(groups)
    n(i) = length(find(grp == groups(i))); 
    
    mean_FR(i)        =  mean(FR_real(find(groups(i) == grp)));
    sem_FR(i)         = (std(FR_real(find(groups(i) == grp))))     / (sqrt(n(i)))   ;
    
    mean_FRsynt(i)    =  mean(FR_synt(find(groups(i) == grp)));
    sem_FRsynt(i)     = (std(FR_synt(find(groups(i) == grp))))     / (sqrt(n(i)))   ;

end

% group order ['CS30   '; 'EpiC   '; 'WTsli32'; 'CS29   '; 'CS30old']

FRs = [FR_real(6:40)';FR_synt(6:40)'];












