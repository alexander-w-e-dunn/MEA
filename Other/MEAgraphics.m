% this script adds or removes custom graphs
% AWE Dunn, Cambridge U., 2020
% INPUTS:
%           add: 0 = removes user-defined options
%                1 = add custom options

function ADgraphics(add)

if add == 0
    
    a = get(groot,'default');
    fields = fieldnames(a);
%     struct2cell(a);
for i = 1 : length(fields)
    set(groot,fields{i,1},'remove')
end
    
elseif add == 1
    
    set(groot,'defaultAxesFontName','Arial')
    set(groot,'defaultAxesFontSize',14)
    set(groot,'defaultAxesLineWidth',1)
    set(groot,'defaultLineLineWidth',1)
%     set(groot,'defaultFigurePosition',[700   435   520   350])
    set(groot,'defaultFigureColor','w'); % white background
    set(groot,'defaultAxesTickDir','out');
    set(groot,'defaultAxesBox','off');
    set(groot,'DefaultAxesTitleFontWeight','Normal')
end

end

% note that some functions override the box and tick dir settings c.f.
% function plotad() customised to override this and imagescad()