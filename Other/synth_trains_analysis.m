% developing a thresholding procedure and testing synth. spike trains ISI
% distriubiton against possion dis.; can then use this to confirm spike
% detection is detecting biological activity; can also use the 95th
% percentile edge weight as the threshold for adjacency matrices for graph
% analysis. Could write a function that opens the adjM as well and saves
% the 95p_threshold in there so it is there as an option.

load('200127_FTDOrg_GrpD_5B_Slice7_mSpikes_3_CTRL_adjM_0.175.mat');

for i = 1:size (adjM2s,3)
    adjMs(:,:,i) = adjM2s(:,:,i) - eye(size(adjM2s(:,:,i)));
end

w = mean(adjMs);
threshold = prctile( reshape(w, [size(w,1), size(w,2) * size(w,3)]) , 95);

real = load('200127_FTDOrg_GrpD_5B_Slice7_mSpikes_3_adjM_0.175.mat');
real_adjM = real.adjM;
bin_adjM = real_adjM;
bin_adjM (find(bin_adjM < threshold)) = 0;
bin_adjM (find(bin_adjM > 0        )) = 1;
bin_adjM = bin_adjM - eye(size(bin_adjM));

% remove NaNs and negative values
bin_adjM (find(isnan(bin_adjM) == 1)) = 0;
bin_adjM (find(bin_adjM < 0        )) = 0;

figure; imagesc(bin_adjM);

nk = sum(sum(bin_adjM)); % num connections
density = nk / ( (size(bin_adjM,1) - 1) * (size(bin_adjM,1) - 2) );
skwnss = skewness(sum(bin_adjM));

for i = 1:size(adjMs,3)
    
    synt_adjM = adjMs (:,:,i);
    bin_adjM2 = synt_adjM;
    bin_adjM2 (find(bin_adjM2 < threshold)) = 0;
    bin_adjM2 (find(bin_adjM2 > 0        )) = 1;
    bin_adjM2 = bin_adjM2 - eye(size(bin_adjM2));
    
    % remove NaNs and negative values
    bin_adjM2 (find(isnan(bin_adjM2) == 1)) = 0;
    bin_adjM2 (find(bin_adjM2 < 0        )) = 0;
    
%     figure; imagesc(bin_adjM2);
    
    nk_v(i) = sum(sum(bin_adjM2)); % num connections
    density_v(i) = nk / ( (size(bin_adjM2,1) - 1) * (size(bin_adjM2,1) - 2) ); %proportion edge density
    
%     % plot degree distributiuon of real(blue) vs synth (red)
%     figure; histogram(sum(bin_adjM)); hold on;  histogram(sum(bin_adjM2));
    
    % more positive means real net was more positively skewed than synth network
    norm_skewness_v(i) = skewness(sum(bin_adjM)) / skewness(sum(bin_adjM2));
    
end

nk_s = mean(nk_v);
density_s = mean(density_v);
norm_skewness = mean(norm_skewness_v);
% proportion of times synth network was less postiively skewed than real
% network:
p_skew = length(find(norm_skewness_v >= 1)) / length(norm_skewness_v);







