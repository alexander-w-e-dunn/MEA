% organoid supp fig
clearvars -except fileName
fileName = '200127_FTDOrg_GrpD_5B_Slice7_mSpikes_3_adjM_0.175.mat';
load(fileName);
adjM(find(channels == 15),:) = 0;
adjM(:,find(channels == 15)) = 0;
% remove negative connections for sttc
adjM(adjM < 0) = 0;
thrOption = 'absolute';
thr = 0.6
badjM1 = threshold_absolute(adjM, thr);
badjM1 = weight_conversion(badjM1, 'autofix');
badjM = weight_conversion(badjM1, 'binarize');
degree = sum(badjM)';
MEAgraphics(1)
badjM = weight_conversion(badjM1, 'binarize');
badjM = badjM + eye(size(badjM));
f1 = figure;
imagesc(badjM)
aesthetics
ylabel('Electrode')
xlabel('Electrode')
set(gca,'TickDir','out');
set(gca, 'FontSize', 14)
fileName1 = fileName
if  contains(fileName1,'_') %remove underscores for title
fileName1(strfind(fileName1,'_'))=' ';
%fileName1=strcat('{',fileName1,'}');
title({strcat(fileName1(1:end-4),' real'),' '});
else
title({strcat(fileName1(1:end-4),' real'),' '});
end
ax = gca;
ax.TitleFontSizeMultiplier = 0.3;
xticks(10:10:60);yticks(10:10:60)
f1.Position = [700   491   329   294];
hold on
for j = 1:length(channels)
plot([j - 0.5 j - 0.5], [0.5 60.5],'LineStyle','-','LineWidth',0.5,'Color',0.4*[1 1 1])
plot([0.5 60.5]       , [j - 0.5 j - 0.5],'LineStyle','-','LineWidth',0.25,'Color',0.4*[1 1 1])
end
hold off
caxis([-1 1])
saveas(f1,strcat(fileName(1:end-4),'_SuppFigAdjM.png'));
saveas(f1,strcat(fileName(1:end-4),'_SuppFigAdjM.epsc'));
close(f1)

% clearvars -except fileName badjM
R = makerandCIJ_und(length(degree),sum(sum(degree))/2);
R = R + eye(size(R));
f2 = figure;
imagesc(R)
aesthetics
ylabel('Electrode')
xlabel('Electrode')
set(gca,'TickDir','out');
set(gca, 'FontSize', 14)
fileName1 = fileName;
if  contains(fileName1,'_') %remove underscores for title
fileName1(strfind(fileName1,'_'))=' ';
%fileName1=strcat('{',fileName1,'}');
title({strcat(fileName1(1:end-4),' surrogate'),' '});
else
title({strcat(fileName1(1:end-4),' surrogate'),' '});
end
ax = gca;
ax.TitleFontSizeMultiplier = 0.3;
f2.Position = [700   491   329   294];
xticks(10:10:60);yticks(10:10:60)
hold on
for j = 1:length(channels)
plot([j - 0.5 j - 0.5], [0.5 60.5],'LineStyle','-','LineWidth',0.5,'Color',0.4*[1 1 1])
plot([0.5 60.5]       , [j - 0.5 j - 0.5],'LineStyle','-','LineWidth',0.25,'Color',0.4*[1 1 1])
end
hold off
colormap(flip(parula))
caxis([-1 1])
saveas(f2,strcat(fileName(1:end-4),'_SuppFigAdjMsurrogate.png'));
saveas(f2,strcat(fileName(1:end-4),'_SuppFigAdjMsurrogate.epsc'));
close(f2)
% caused low level graphics error so i used command:
% opengl('save', 'software')
% had received this:
% Warning: MATLAB previously crashed due to a low-level graphics error. To
% prevent another crash in this session, MATLAB is using software OpenGL instead
% of using your graphics hardware. To save this setting for future sessions, use
% the opengl('save', 'software') command. For more information, see Resolving
% Low-Level Graphics Issues. 
% > In matlab.graphics.internal.initialize (line 15) 

R = R - eye(size(R));
badjM = badjM - eye(size(R));
comp = badjM - R;
f3 = figure;
colormap(parula)
imagesc(comp)
aesthetics
ylabel('Electrode')
xlabel('Electrode')
set(gca,'TickDir','out');
set(gca, 'FontSize', 14)
fileName1 = fileName;
if  contains(fileName1,'_') %remove underscores for title
fileName1(strfind(fileName1,'_'))=' ';
%fileName1=strcat('{',fileName1,'}');
title({strcat(fileName1(1:end-4),' real-surrogate'),' '});
else
title({strcat(fileName1(1:end-4),' real-surrogate'),' '});
end
ax = gca;
ax.TitleFontSizeMultiplier = 0.3;
f3.Position = [700   491   329   294];
xticks(10:10:60);yticks(10:10:60)
hold on
for j = 1:length(channels)
plot([j - 0.5 j - 0.5], [0.5 60.5],'LineStyle','-','LineWidth',0.5,'Color',0.4*[1 1 1])
plot([0.5 60.5]       , [j - 0.5 j - 0.5],'LineStyle','-','LineWidth',0.25,'Color',0.4*[1 1 1])
end
hold off
saveas(f3,strcat(fileName(1:end-4),'_SuppFigAdjMreal-surrogate.png'));
saveas(f3,strcat(fileName(1:end-4),'_SuppFigAdjMreal-surrogate.epsc'));
close(f3)

surrdegree = sum(R);
f4 = figure;
f4.Position = [700   435   422   350];
c = parula; yell = c(end,:); blau = c(1,:);
histogram(degree,'FaceColor',yell,'LineWidth',2);
aesthetics
ylim([0 15]); 
ys = ylim; xs = xlim; xts = xticks; yts = yticks;
xlabel('node degree');
ylabel('frequency');
legend('real','box','off')
saveas(f4,strcat(fileName(1:end-4),'_SuppFighistReal.png'));
saveas(f4,strcat(fileName(1:end-4),'_SuppFighistReal.epsc'));
close(f4)

surrdegree = sum(R);
f5 = figure;
f5.Position = [700   435   422   350];
histogram(surrdegree,'FaceColor','b','LineWidth',2);
aesthetics
ylim(ys); 
xlim(xs); xticks(xts); yticks(yts);
xlabel('node degree');
ylabel('frequency');
% legend('surrogate','box','off')
legend('synth.','box','off')
saveas(f5,strcat(fileName(1:end-4),'_SuppFighistSurr.png'));
saveas(f5,strcat(fileName(1:end-4),'_SuppFighistSurr.epsc'));
close(f5)

