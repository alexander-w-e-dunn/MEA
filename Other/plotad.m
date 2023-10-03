 function varargout = plotad(varargin)
    % Override builtin plot and honors the following settings of the
    % graphics root that builtin plot for some reason ignores:
    %   'DefaultAxesBox'
    %   'DefaultAxesTickDir'
    h=builtin('plot',varargin{:});
    set(get(h,'Parent'),'Box',get(groot,'DefaultAxesBox'));
    set(get(h,'Parent'),'TickDir',get(groot,'DefaultAxesTickDir'));
    varargout{1}=h;  
end