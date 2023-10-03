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

culture = 1;
MEAenergy=HDMEA1449energy;

for rule = 2 : size(MEAenergy,2)
    
    % get lowest energy value for each rule
%     optimum_params = find 
    
    nColumns = 3;
    subplot(ceil(size(MEAenergy,2)/nColumns) , nColumns, rule-1)
    imagesc(reshape(MEAenergy(culture,rule,:),sqrt(size(MEAenergy,3)) , sqrt(size(MEAenergy,3)) ))
    set(gca, 'XTick', [1 sqrt(size(MEAenergy,3))], 'XTickLabel', [-7 7])
    set(gca, 'YTick', [1 sqrt(size(MEAenergy,3))], 'YTickLabel', [-7 7])
    ylabel('\eta')
    xlabel('\gamma')
    title(modeltype(rule))  
    caxis([0 1])
end

f1.Position = [680   137   560   841];

% imagesc(reshape(HDMEA1449energy(1,3,:),10,10));
% xticklabels