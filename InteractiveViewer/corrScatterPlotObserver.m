% obs = corrScatterPlotObserver(desMtx, wh_of_interest, onsets)
%
% Example:
% desMtx = trialmodels.ratings;
% wh_of_interest = 1;
% images = spm_get(Inf,'*img','Select images for timeseries plot');
% % < Select trial-level magnitude images for a subject >
%
%
% iv = InteractiveViewer('UseExistingGraphicsWindow', 1, 'AssociatedVols', images, 'IVObserver', corrScatterPlotObserver(desMtx, 1));
% set(gca,'FontSize',16)
% xlabel('Pain ratings (centered)');
% ylabel('fMRI Response magnitude');
% %
%
% One col. of ratings, nothing to adjust
% iv = InteractiveViewer('UseExistingGraphicsWindow', 1, 'AssociatedVols', ...
% images, 'IVObserver', corrScatterPlotObserver(desMtx(:,1), 1));

function obs = corrScatterPlotObserver(desMtx, wh_of_interest)

    % Do things here that are not changed every time you click
    % Any variable I define here will be avail. in the inline, and thus in
    % the interactive viewer.
    
% %     [n, k] = size(desMtx);
    
% %     px = pinv(desMtx);
    
    desMtx_noint = desMtx;
    desMtx_noint(:, wh_of_interest) = [];
    
    predictor = desMtx(:, wh_of_interest);
    
    px_noint = pinv(desMtx_noint);
    
    
    obs = IVObserver(@plot_);

    % This inline function is run whenever you click
    function plot_(iv, mmPos, voxPos)
        ts = get(iv, 'CurrentVolsTs');
        
% %         ndata = length(ts);
% %         if ndata < n
% %             disp('Warning! Timeseries has fewer data points than model!')
% %             desMtx = desMtx(1:ndata, :);
% %             px = px(:, 1:ndata);
% %             
% %             desMtx_noint = desMtx_noint(1:ndata, :);
% %             px_noint = px_noint(:, 1:ndata);
% %         end

        if ~iscol(ts), ts = ts'; end
% %         b = px * ts;
        
        % predictor values, adjusted for effects of no interest
        adj_pred = predictor - desMtx_noint * px_noint * predictor;
        
        % adjusted data, removing vars of no interest
        ts_adj = ts - desMtx_noint * px_noint * ts;
        
        hold off;
        plot(adj_pred, ts_adj, 'bo', 'MarkerFaceColor', [.5 .5 .5]);
        
        refline
        
        
        title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f', mmPos(1), mmPos(2), mmPos(3)), 'FontSize', 14);
        
        % Return data to workspace
        disp('Assigning braindata with image values in workspace.');
        assignin('base', 'braindata', ts);
        
    end
end

