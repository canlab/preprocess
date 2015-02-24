function [f1,colors,axh] = tor_fig(varargin)
% function [f1,colors,axh] = tor_fig(varargin)
%
% Creates a figure, hold on, fontsize 18
% 
% Var args are number of rows, then columns for subplot

f1 = figure('Color','w'); set(gca,'FontSize',18),hold on
colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};

if length(varargin) > 0
    i = max(1, varargin{1});
    j = max(1, varargin{2});
    np = max(1, i * j);
    
    for k = 1:np
        axh(k) = subplot(i,j,k); 
        set(gca,'FontSize',18),hold on
    end
    axes(axh(1));
end
    
return