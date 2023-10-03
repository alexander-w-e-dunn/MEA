function aesthetics() 
set(gca, 'box', 'off') % remove borders
set(gcf,'color','w'); % white background
set(gca, 'TickDir', 'out')
set(gca,'FontSize',14);
set(gca,'FontName','Arial');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'LineWidth',2);
ax = gca;
ax.TitleFontWeight = 'Normal';
try
    ax.Colorbar.TickDirection = 'out';
    ax.Colorbar.FontName = 'Arial';
    ax.Colorbar.FontSize = 14;
    ax.Colorbar.Label.FontSize = 16;
    ax.Colorbar.Box = 'off';
catch
end 

try
    ax.Legend.FontName = 'Arial';
    ax.Legend.FontSize = 14;
    ax.Legend.FontWeight = 'normal';
    ax.Legend.Title.FontWeight = 'normal';
    ax.Legend.Box = 'off';
catch
end
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
%     'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:500:2500, ...
%     'LineWidth', 1)

% set(gca, 'Box', 'off', 'TickDir', 'out',  ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on')
set(gca, 'YGrid','off','XGrid','off')
end 