files = dir('MPT200107_1A_DIV28*0.0627*adjM_0.05.mat');
load(files(1).name);
adjM1 = adjM - 6 * eye(size(adjM));
adjM1(find(isnan(adjM1))) = -5;
f1=figure; imagesc(adjM1);
w = find(adjM1 >-5);
thr = prctile(adjM1(w),60);
adjM2 = adjM1;
adjM2(find(adjM2 < thr)) = 0; 
adjM3 = adjM2;
adjM2(find(adjM2 >= thr)) = 1;
f2=figure; imagesc(adjM2);
f3 = figure; histogram(adjM3(find(adjM3 >= thr)));
f4 = figure; histogram(adjM1(w));