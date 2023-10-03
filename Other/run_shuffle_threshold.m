clear all

filename = '200127_FTDOrg_GrpA_3A_Slice8_mSpikes_3.mat';

cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
load(filename)
events = full(mSpikes);
% events = full(mSpikes)';
load(strcat(filename(1:end-4),'_adjM.mat'));
% repNum = 100;
repNum = 100;
randMethod = 'circular'
adjMmethod = 'tileCoef'
fs = 25000
rec_s = length(events)/fs
downSample = rec_s*1000
ds_rate = downSample/length(events)
sync_win_s = 0.175;
lag = ds_rate * sync_win_s;
addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output\RF_scripts'))
% tic 
[CI,adjMm,adjMci,shuffled_spikeMat_sparse_all,adjMi] = AdjMProbabilisticThreshMEA(events,adjM,repNum,randMethod,adjMmethod,downSample,lag);
% toc

% figure;subplot(1,2,1);imagesc(adjM);subplot(1,2,2);imagesc(adjMci)
figure;
subplot(1,2,1);imagesc(CI(:,:,1));colorbar; title('CI Lower Bounds'); caxis([min(CI(:)) , max(CI(:))])
subplot(1,2,2);imagesc(CI(:,:,2));colorbar; title('CI Upper Bounds'); caxis([min(CI(:)) , max(CI(:))])
% aesthetics
% clearvars -except events adjM repNum randMethod adjMmethod downSample lag
figure; 
subplot(1,2,1); imagesc(mean(adjMi,3));colorbar; caxis([0 1]); title('Mean Synthetic')

threshold = 95;

for row = 1 : size(adjMi,1)
    for column = 1 : size(adjMi,2)
        if row ~= column && row ~= 15 && column ~= 15
            rand_edges = squeeze(adjMi(row,column,:));
            thr(row,column) = prctile(rand_edges,threshold);
            msyn(row,column) = mean(rand_edges); 
        end
    end
end

subplot(1,2,2); imagesc(thr);colorbar; caxis([0 1]); title('95^t^h Percentile Synthetic')

save(strcat(filename(1:end-4),'_synth_data2.mat'),'adjMi','shuffled_spikeMat_sparse_all','filename','-v7.3')

%%
figure; n = 3;
for i = 1 : n
    subplot(1,n,i)
    imagesc(adjMi(:,:,i))
    title({['Randomised adj.m. #', num2str(i)]})
end

%% degree
ba = zeros(60);
ba(find(adjMci>0)) = 1;
% ba(isnan(ba)) = 0;
% ba = weight_conversion(ba,'binarize');
ba = weight_conversion(ba,'autofix');
figure;subplot(1,2,1);imagesc(triu(ba));
k = sum(triu(ba));
subplot(1,2,2);histfit(k,20,'kernel')