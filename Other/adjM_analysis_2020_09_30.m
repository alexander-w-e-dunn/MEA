
%% load adjM, get vector of weights, plot histogram
clear all
close all
MEAgraphics(1)
addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output\2019_03_03_BCT'))
load('C:\Users\alexd\Dropbox (Cambridge University)\NOG MEA Data\MEA Data Mecp2 Project Jan 2019-\Generative_model_format\MEA_matrices_EpiC_human_iPSC_all_ages.mat')
a=adjM_all(:,:,2);
a=weight_conversion(a,'autofix');
f1 = figure; f1.Position  = [ 392   208   839   623];
subplot(2,2,1)
imagesc(a+eye(size(a)));aesthetics;cb=colorbar;cb.Visible='off';
ytsl = yticklabels; yts = yticks;
xticks(yts); xticklabels(ytsl);
title('(a)', 'FontSize', 20, 'HorizontalAlignment', 'left','position',[-18 0]);
% nansum(a)
w=triu(a);
% imagesc(w);
b=ones(60);
b=b.*10;
b=triu(b,1);
% subplot(2,3,2);
% imagesc(b)
c=w+b;
subplot(2,2,2)
imagesc(c-b); yticklabels([]);xticks(yts); xticklabels(ytsl);aesthetics
cb=colorbar; aesthetics
title('(b)', 'FontSize', 20, 'HorizontalAlignment', 'left','position',[-18 0]);

w2=c(:);
w2;
w2=w2(w2~=0);
w2=w2-10;
subplot(2,2,[3 4])
h=histfit(w2,20,'kernel'); aesthetics; xlabel('Correlation'); ylabel ('Frequency');
h(1).EdgeColor = [1 1 1]; h(1).FaceColor = 0*[1 1 1];   
h(2);
ylims = ylim;
violinx = h(2).YData;
violiny = h(2).XData;
violinylim = round(min(w2),1) : 0.1 : round(max(w2),1) ;
hold on; findpeaks(violinx,violiny); aesthetics; xlabel('Correlation'); ylabel ('Frequency'); ylim(ylims);
[pks, pklocs] = findpeaks(violinx,violiny);
cb=colorbar;cb.Visible='off';
title('(c)', 'FontSize', 20, 'HorizontalAlignment', 'left','position',[-0.8 310]);

%% get ctrl adjM
% if ctrl adjM file doesn't already exist, call function to get 100 adjMs


%% find degree cut off
% find trough immediately before the peak
% could use findpeaks on negative vector (findpeaks(-violinx,violiny)) and
% use the first trough or could look at the curve between the first 2 peaks
% and find the minimum value
[trghs, trghlocs] = findpeaks(-violinx,violiny);

% determine cut off point for high degree nodes
if pks > 1
    cutoff = abs(trghs(1));
    %add variable for high degree nodes highdegnodes = find nodes > cutoff
    % plot cutoff
    hold on
    plot(trghlocs(1),cutoff,'or','MarkerSize',10)
end

