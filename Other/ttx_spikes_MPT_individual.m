clear all; close all

filesnames = {'MPT190403_6B_DIV50TTXBEFORE_cSpikes_L-0.0627_RF1.mat';... % this will be on y axis
    'MPT190403_6B_DIV50TTXBEFORE_mSpikes_4.5_RF1.mat'};                  % this will be on x axis
% fileName = 'MPT190403_6B_DIV50TTXBEFORE_cSpikes_L-0.0627_RF1.mat';
% fileName = 'MPT190403_6B_DIV50TTXBEFORE_mSpikes_4.5_RF1.mat';
% filesnames = {'191210_slice1_DIV_g04_2018_2min_cSpikes_L-0.25071_RF1.mat';... % this will be on y axis
%     '191210_slice1_DIV_g04_2018_2min_mSpikes_3_RF1.mat'};                  % this will be on x axis
% postfiles = {'191210_slice1_DIV_g04_2018_TTX_2min.mat';... % this will be on y axis
%     '191210_slice1_DIV_g04_2018_TTX_2min.mat'};

for file = 1:length(filesnames)
    fileName = filesnames{file};
    load(fileName)
    mean(thresholds)
    try
        pre = full(sum(mSpikes));
        preSparseM = mSpikes;
        pre_length = length(mSpikes);
    catch
        pre = full(sum(cSpikes));
        pre_length = length(cSpikes);
        preSparseC = cSpikes;
    end
    % load( [ fileName(1:strfind(lower(fileName),'ttx')-1) , fileName(strfind(lower(fileName),'spikes')-2 : end) ] )
    if ~isempty(strfind(fileName,'MPT'))
        postfile = [ fileName(1:strfind(lower(fileName),'ttx')-1) , 'TTXAFTER0001.mat' ];
    elseif ~isempty(strfind(fileName,'slice'))
        postfile = postfiles{file};
    end
    L=0; refPeriod_ms = num2str(fileName(strfind(fileName,'RF')+2 : end-4));
    
    spikeparams = fileName( strfind(fileName,'Spikes') : end-4);
    % refPeriod_ms = str2num( spikeparams (strfind(spikeparams,'RF')+2 : end) );
    if contains(fileName,'cSpikes')
        method = 'cwt';
        L = str2double( spikeparams(strfind(spikeparams,'_L')+2 : strfind(spikeparams,'_RF')-1) );
    elseif contains(fileName,'mSpikes')
        method = 'Manuel';
    end
    
    if ~~exist('mSpikes')
        
        load( postfile );
        % remove grounded
        thresholds(find(thresholds > -9)) = -900; %original value was -10
%         thresholds(find(thresholds >= thresholds(find(channels==15)))) = -900; %use ref channel to find grounded
        progressbar('elecs')
        for elec = 1:length(thresholds)
            
            multiplier = thresholds(elec);
            [spikeMatrix(:, elec), finalData(:, elec), thresholds(elec)] = detectSpikes(dat(:, elec), 'abs', multiplier,L,refPeriod_ms);
            progressbar(elec/length(thresholds))
        end
        plotcol = 'r';
        trajcol = [1 0 0 0];
        postSparseM = sparse(spikeMatrix);
        
    elseif ~~exist('cSpikes')
        
        postdata = load( strcat( postfile(1:end-4) , '_cSpikes_L' , ...
            num2str(L) , '_RF', num2str(refPeriod_ms) , '.mat' ) );
        spikeMatrix = full(postdata.cSpikes);
        postSparseC = sparse(spikeMatrix);
        plotcol = 'b';
        trajcol = [0 0 1 0];
    end
    
    recTimeRatio = pre_length / length(spikeMatrix) ;
    post = round( sum(spikeMatrix) * recTimeRatio);
%     a = pre-post
%     a(find(channels == 22))
%     pre(find(channels == 22))
%     post(find(channels == 22))
% remove channels that are grounded in pre ttx rec
    post(pre<=1) = 0;
%     mean(pre-post)
%     % histogram(pre-post)
%     mean(pre)
%     mean(post)
%     mean(pre) - mean(post)
    
    f1 = figure;
    jitval = 0.1;
    y = log10( [pre' , post'] ); y(find(y==-Inf)) = 0;
    % x = [ones(size(pre))' , 2* ones(size(post))'];
    x = [1 + 1.2*jitval, 2 - 1.2*jitval];
    plot(x , y , '--r','Visible','off'); hold on
    xlimits = [0.5 2.5];
    xlim(xlimits)
    xticks([1 2]); xticklabels({'baseline';'TTX'});
        % move axis so you can see 0 values...
    ylim([-0.2 4])
    yticks(log10([1 5 10 50 100 500 1000 5000 10000]));
%     ylim([-0.2 5])
%     yticks(log10([1 5 10 50 100 500 1000 5000 10000 50000 100000]));
    yticklabels(10.^(yticks))
    ylabel('spike count')
    % plines = plot(xlimits,[yticks ; yticks],'-','Color',[0 0 0 .2],'LineWidth',2);
    gridlines = plot([0.52 2.48],[yticks ; yticks],'-','Color',0.75*[1 1 1],'LineWidth',1.5);
    trajectories = plot(x , y , ':','LineWidth',1.5,'Color',trajcol + [0 0 0 .2]);
    SEMs = std(y) ./ sqrt(length(y));
    x1 = x; y1 = y;
    x = [ones(size(pre)) , 2* ones(size(post))];
    y = reshape(y,[1 numel(y)]);
    ind_points = scatter(x,y,'jitter', 'on', 'jitterAmount',jitval,'MarkerEdgeColor',plotcol,'LineWidth',1.5,...
        'MarkerEdgeAlpha',.5);
    aesthetics
    bx = gca;
    bx.TickDir = 'out';
    bx.TickLength = 0.02 * [1 1];
    bx.LineWidth = 2;
    bx.FontName = 'Arial';
    bx.FontSize = 18;
    f1.Position = [1030 270 424 420];
    % drawaxis(bx, 'x', -0.2, 'movelabel', 1) % download from https://uk.mathworks.com/matlabcentral/answers/uploaded_files/2281/drawaxis.m
    ebars = errorbar([1 2],mean(y1),SEMs,'LineStyle','none','Color','k','LineWidth',2,'CapSize',18);
    e_mean_squares = plot([1 2],mean(y1),'s','LineStyle','none','Color','k','LineWidth',2,'MarkerFaceColor','k')
    
    fileName1 = fileName;
    if  contains(fileName,'_') %remove underscores for title
        fileName1(strfind(fileName1,'_'))=' ';
        %fileName1=strcat('{',fileName1,'}');
        title({fileName1(1:end-4),' TTX experiment',' ',' '});
    else
        title({fileName1(1:end-4),' TTX experiment',' ',' '});
    end
    %make title font smaller
    bx.TitleFontSizeMultiplier = 0.45;
    
    %save raster as PNG
    saveas(f1,strcat(fileName(1:end-4),'_TTX_experiment',lower(method(1)),spikeparams,'.png'));
    if ~~exist('f1')
        close(f1);
    end
    
    clear mSpikes cSpikes
    
end
clearvars -except preSparseM preSparseC postSparseM postSparseC channels ...
    thresholds recTimeRatio filesnames method spikeparams

% corrected_channels = find(thresholds == -900);

befvals = [full(sum(preSparseC))  ;  full(sum(preSparseM))];
aftvals = [full(sum(postSparseC)) ;  full(sum(postSparseM))] .* recTimeRatio;

% % correction so that when plotting log of spikes, 0s will be plotted as
% % they are currently excluded due to being -Inf
befvals(find(befvals == 0)) = 1;
aftvals(find(aftvals == 0)) = 1;
% befvals(:,corrected_channels) = 0;
% aftvals(:,corrected_channels) = 0;

figBEF = figure;
figBEF.Position = [680 558 1120 420]
subplot(1,2,1)
markersize = 50;
upperlim = log10(10000);
% upperlim = log10(100000);
xlim([0 upperlim]) ; ylim([0 upperlim]); hold on
plot([0.1 upperlim-0.1],[0.1 upperlim-0.1],':','Color',0.7*[1 1 1],'LineWidth',3);
p2 = scatter(log10(befvals(2,:)) , log10(befvals(1,:)),markersize, 'LineWidth',2,'MarkerEdgeColor','k');
xticks(0:upperlim);
xticklabels(10.^xticks);
yticks(0:upperlim);
yticklabels(10.^yticks);
aesthetics
bx = gca;
bx.TickDir = 'out';
bx.TickLength = 0.02 * [1 1];
bx.LineWidth = 2;
bx.FontName = 'Arial';
bx.FontSize = 18;
xlabel({'spike count','(threshold method)'});
ylabel({'spike count','(template method)'});

% thr1 = load(filesnames{2},'thresholds');thresholds1 = thr1.thresholds; thresholds1(find(thresholds1 >= thresholds1(find(channels==15)))) = -900;
% inactive_vec = find(thresholds1 == -900); %i.e. a low ampl. vec.
inactive_vec = find(befvals(1,:) < 9);
p3 = scatter(log10(befvals(2,inactive_vec)) , log10(befvals(1,inactive_vec)),...
    markersize,'LineWidth',2,'MarkerEdgeColor','m');
legendinho = legend([p2 p3],{'active','inactive'},'Location','southeast');
title(legendinho,'status','FontSize',12)
% xtickangle(bx,45); legend off;
% M = findobj(legendinho,'type','patch');
% set(M, 'Markersize', markersize); %// set value as desired

fileName = filesnames{1,1};
fileName1 = fileName;
fileName2 = filesnames{2,1};
if  contains(fileName,'_') %remove underscores for title
    fileName1(strfind(fileName1,'_'))=' ';
    fileName2(strfind(fileName2,'_'))=' ';
    %fileName1=strcat('{',fileName1,'}');
    title({fileName1(1:end-4),fileName2(1:end-4),' TTX experiment: BEFORE TTX',' ',' '});
else
    title({fileName1(1:end-4),fileName2(1:end-4),' TTX experiment: BEFORE TTX',' ',' '});
end
%make title font smaller
bx.TitleFontSizeMultiplier = 0.45;




subplot(1,2,2)
markersize = 50;
% upperlim = log10(10000);
xlim([0 upperlim]) ; ylim([0 upperlim]); hold on
plot([0.1 upperlim-0.1],[0.1 upperlim-0.1],':','Color',0.7*[1 1 1],'LineWidth',3);
p4 = scatter(log10(aftvals(2,:)) , log10(aftvals(1,:)),markersize, 'LineWidth',2,'MarkerEdgeColor','k');
xticks(0:upperlim);
xticklabels(10.^xticks);
yticks(0:upperlim);
yticklabels(10.^yticks);
aesthetics
bx = gca;
bx.TickDir = 'out';
bx.TickLength = 0.02 * [1 1];
bx.LineWidth = 2;
bx.FontName = 'Arial';
bx.FontSize = 18;
xlabel({'spike count','(threshold method)'});
ylabel({'spike count','(template method)'});

% inactive_vec = find(thresholds == -900); %i.e. a low ampl. vec.
inactive_vec = find(befvals(1,:) < 9);
p5 = scatter(log10(aftvals(2,inactive_vec)) , log10(aftvals(1,inactive_vec)),...
    markersize,'LineWidth',2,'MarkerEdgeColor','m');
legendinho1 = legend([p4 p5],{'active','inactive'},'Location','southeast');
title(legendinho1,'pre-TTX status','FontSize',12)
% xtickangle(bx,45); legend off;
% M = findobj(legendinho,'type','patch');
% set(M, 'Markersize', markersize); %// set value as desired

fileName = filesnames{1,1};
fileName1 = fileName;
fileName2 = filesnames{2,1};
if  contains(fileName,'_') %remove underscores for title
    fileName1(strfind(fileName1,'_'))=' ';
    fileName2(strfind(fileName2,'_'))=' ';
    %fileName1=strcat('{',fileName1,'}');
    title({fileName1(1:end-4),fileName2(1:end-4),' TTX experiment: AFTER TTX',' ',' '});
else
    title({fileName1(1:end-4),fileName2(1:end-4),' TTX experiment: AFTER TTX',' ',' '});
end
%make title font smaller
bx.TitleFontSizeMultiplier = 0.45;



%save comparison as PNG
fileName2 = filesnames{2,1};
fileName1 = filesnames{1,1};
spikeparams1 = fileName1( strfind(fileName1,'Spikes') : end-4);
spikeparams2 = fileName2( strfind(fileName2,'Spikes') : end-4);

%% SAVE 
saveas(figBEF,strcat(fileName1(1:strfind(fileName1,'Spikes')-3),'_',...
    fileName1(strfind(fileName1,'Spikes')-1 : end-4),'_and_',...
    fileName2(strfind(fileName2,'Spikes')-1 : end-4),...
    '_TTX_comparison.png'));
if ~~exist('figBEF')
    close(figBEF);
end


