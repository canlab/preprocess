function obs = tsAncovaPlotObserver(plotAgainst, controlFor, groupBy, varargin)
    % obs = tsAncovaPlotObserver(plotAgainst, controlFor, groupBy, varargin)
    
    if length(varargin) > 0
        cl = varargin{1};
    end


    X = controlFor;
    n = size(X,1);
    iV = eye(n);        % would be whitening matrix, but leave as identity
    
    R = (eye(n) - X * inv(X'*iV*X)*X'*iV);        % Residual inducing matrix
       
    groupValues = unique(groupBy);
    
    colors = {[1 1 .3] [.3 .3 1]};
    
    % create a function handle with ?
    obs = IVObserver(@ts_plot_);
    

% This inline function is run whenever you click
% It has access to all vars in main script, which become embedded in the
% function handle

function ts_plot_(iv, mmPos, voxPos)
    ts = get(iv, 'CurrentVolsTs');
    
    % OLS
    ts_adj = R * ts;
    
    
    % robust
    [b, stat] = robustfit([plotAgainst groupBy controlFor], ts);
    
    cla
    hold on;
    for i = 1:length(groupValues)
        wh = groupBy == groupValues(i);
        plot(plotAgainst(wh), ts_adj(wh), 'o', 'Color', [.2 .2 .2], 'MarkerFaceColor', colors{i});
        refline
        
    end
    
    title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f', mmPos(1), mmPos(2), mmPos(3)), 'FontSize', 14);
    
    
    % text output
    % ------------------------------------------------------
        
    fprintf(1, 'Robust output for this voxel:\n-------------------------------------\n')
    fprintf('Parameter\t\test.\tSE\tt(%3.0f)\tp\n', stat.dfe)

    i = 1; fprintf(1, 'Intercept\t\t%3.2f\t%3.2f\t%3.2f\t%3.4f\t\n', b(i), stat.se(i), stat.t(i), stat.p(i) )
    i = 2; fprintf(1, 'plotAgainst\t\t%3.2f\t%3.2f\t%3.2f\t%3.4f\t\n', b(i), stat.se(i), stat.t(i), stat.p(i) )
    i = 3; fprintf(1, 'groupingVar\t\t%3.2f\t%3.2f\t%3.2f\t%3.4f\t\n', b(i), stat.se(i), stat.t(i), stat.p(i) )
    
    for i = 4:size(controlFor, 2) + 3
        fprintf(1, 'Covariate\t\t%3.2f\t%3.2f\t%3.2f\t%3.4f\t\n', b(i), stat.se(i), stat.t(i), stat.p(i) )
    end
    fprintf(1, '\n')
    
    
    
    % cluster stuff, if entered
    % ------------------------------------------------------
    
    if exist('cl','var')
        [clusters, wh] = find_closest_cluster(cl, mmPos');
        fprintf('Closest cluster is: %3.0f, Voxels: %3.0f, [x y z] = %3.0f %3.0f %3.0f\n', wh, clusters.numVox, mmPos(1), mmPos(2), mmPos(3))
    end
    
end


end