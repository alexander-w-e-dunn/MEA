

%% set up
clear all; close all
% set directories and plotting graphics
data_directory = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'))
MEAgraphics(0); MEAgraphics(1)
%% load data
cd(data_directory);
files = dir('*slice*.mat');                  
files = files(~~contains({files.name}, 'CTRL'));
% load excel file containing a list of slices and lines
groupdata = readcell('FTD_organoid_groups');
% find the first EpiC slice
Index = find(contains(groupdata(:,3),'EpiC'));
% load(strcat(groupdata{Index(1)},'_mSpikes_3_CTRL_spikeMat.mat'));
load(strcat(groupdata{Index(1)},'_mSpikes_3_CTRL_adjM_0.175.mat'));
load(strcat(groupdata{Index(1)},'_mSpikes_3_adjM_0.175.mat'));
%% get threshold for each edge
threshold = 80;

for row = 1 : size(adjM2s,1)
    for column = 1 : size(adjM2s,2)
        if row ~= column && row ~= 15 && column ~= 15
            rand_edges = squeeze(adjM2s(row,column,:));
            thr(row,column) = prctile(rand_edges,threshold);
            msyn(row,column) = mean(rand_edges); 
        end
    end
end

% remove edges less than OR EQUAL TO threshold i.e. if threshold is 1, edge
% will be removed; 
thr_adjm = adjM; thr_adjm(adjM<=thr)=0; 
% thr_adjm = adjM; thr_adjm(adjM<thr)=0; % this option removes edges <
% threshold i.e. keeps edges with a threshold of 1
f1 = figure; 
subplot(3,2,1)
imagesc(adjM);      caxis([-0 1]); colorbar; aesthetics; title('Observed network');
subplot(3,2,4)
imagesc(thr);       caxis([-0 1]); colorbar; aesthetics; title(strcat(num2str(threshold),'th percentile'));
subplot(3,2,3)
imagesc(msyn);      caxis([-0 1]); colorbar; aesthetics; title('Mean randomisation');
subplot(3,2,5)
imagesc(thr_adjm);  caxis([-0 1]); colorbar; aesthetics; title('Thresholded matrix');
subplot(3,2,6)
imagesc(adjM>thr);  caxis([-0 1]); colorbar; aesthetics; title('Binarised matrix');
% use below to include edges with a threshold of 1
% imagesc(adjM>=thr);  caxis([-0 1]); colorbar; aesthetics; title('Binarised matrix');
f1.Position = [2          42        1002         954]
% thr_adjm = adjM; 
% thr_adjm(thr_adjm < thr

% plot distribution across all edges and iterations
idmat = reshape( repmat( eye(size(adjM)),[1 1 100]), [1, size(adjM,1) * size(adjM,2) * size(adjM2s,3)] ) ;
allrandedges = reshape(adjM2s , [1, size(adjM2s,1) * size(adjM2s,2) * size(adjM2s,3)] );
allrandedges(find(idmat==1)) = NaN;
% allrandedges(allrandedges == 1) = NaN;
subplot(3,2,2)
histfit(allrandedges,50,'kernel')
aesthetics
xlabel('Correlation (STTC)')
ylabel('Frequency')
title('Randomised edges')