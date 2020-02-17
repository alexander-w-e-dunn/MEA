% re order an adjacency matrix
%{
Take an MEA adjM (60 x 60) and reorder according to the channels variable.
Put in order 12, 13, 14 etc. (Channels IDs)
Or may need to change to 21,31,41, etc. 
%}
clear all
files = dir('*191209_FTD_slice7_*cSpikes_L0_adjM*.mat');  % make sure all the spike matrix files in the 
files = files(~contains({files.name}, 'reord'));%remove unwanted files
files = files(~contains({files.name}, 'm_cS'));%remove unwanted files

progressbar
for i=1:length(files)
    tic
    load(files(i).name)
    %{
    topo_order_vector=[find(channels==12),find(channels==13),find(channels==14),...
        find(channels==15),find(channels==16),find(channels==17),find(channels==21),...
        find(channels==22),find(channels==23),find(channels==24),find(channels==25),...
        find(channels==26),find(channels==27),find(channels==28),find(channels==31),...
        find(channels==32),find(channels==33),find(channels==34),find(channels==35),...
        find(channels==36),find(channels==37),find(channels==38),find(channels==41),...
        find(channels==42),find(channels==43),find(channels==44),find(channels==45),...
        find(channels==46),find(channels==47),find(channels==48),find(channels==51),...
        find(channels==52),find(channels==53),find(channels==54),find(channels==55),...
        find(channels==56),find(channels==57),find(channels==58),find(channels==61),...
        find(channels==62),find(channels==63),find(channels==64),find(channels==65),...
        find(channels==66),find(channels==67),find(channels==68),find(channels==71),...
        find(channels==72),find(channels==73),find(channels==74),find(channels==75),...
        find(channels==76),find(channels==77),find(channels==78),...
        find(channels==82),find(channels==83),find(channels==84),find(channels==85),...
        find(channels==86),find(channels==87)];
    above is ID in numerical order
    R script coord variable is in below order
    %}
        topo_order_vector=[find(channels==17),find(channels==16),find(channels==15),...
        find(channels==14),find(channels==13),find(channels==12),find(channels==28),...
        find(channels==27),find(channels==26),find(channels==25),find(channels==24),...
        find(channels==23),find(channels==22),find(channels==21),find(channels==38),...
        find(channels==37),find(channels==36),find(channels==35),find(channels==34),...
        find(channels==33),find(channels==32),find(channels==31),find(channels==48),...
        find(channels==47),find(channels==46),find(channels==45),find(channels==44),...
        find(channels==43),find(channels==42),find(channels==41),find(channels==58),...
        find(channels==57),find(channels==56),find(channels==55),find(channels==54),...
        find(channels==53),find(channels==52),find(channels==51),find(channels==68),...
        find(channels==67),find(channels==66),find(channels==65),find(channels==64),...
        find(channels==63),find(channels==62),find(channels==61),find(channels==78),...
        find(channels==77),find(channels==76),find(channels==75),...
        find(channels==74),find(channels==73),find(channels==72),find(channels==71),...
        find(channels==87),find(channels==86),find(channels==85),find(channels==84),...
        find(channels==83),find(channels==82)];
    
    %below fixes adjM in R script but causes self connections in network
    %plot
    
    %adjM=adjM(flip(topo_order_vector),topo_order_vector);
    %NOTE the above line was tested at one stage because the adjacency
    %matrix in R was plotted in this order, however, the network plot still
    %requires the matrix to be in the below order thus the adjacency matrix
    %figure in R is incorrect, but the network plot will be correct if you
    %use the below order. If you use the above order, the R adjacency
    %matrix plot will be correct but the network plot will be incorrect    
    adjM=adjM(topo_order_vector,topo_order_vector);
    new_channels=channels(topo_order_vector);


    
    
    %figure
    %imagesc(adjM);
    %plotAdj(adjM,1:60);
    %adjM(find(isnan(adjM)))=0;
    %sum(sum(adjM));
    fileName = strcat(files(i).name(1:end-4), '_reord.mat'); 
        % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
        % save(fileName, 'adjM','channels','new_channels');
        save(fileName, 'adjM','channels','new_channels');
    %toc;
        progressbar(i/length(files));
    clearvars -except i files
end

%% plot adjM as it will look in R
load('191209_FTD_slice7_GroupC_1_cSpikes_L0_adjM_reord.mat')
adjM=adjM(flip(1:length(adjM)),1:length(adjM)');
adjM=adjM-flip(eye(size(adjM)));
figure;imagesc(adjM);colorbar
yticks([0 10 20 30 40 50 60])
yticklabels(flip(yticklabels))
yticks([1 10 20 30 40 50 60])