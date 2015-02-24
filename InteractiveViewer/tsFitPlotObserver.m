% obs = tsFitPlotObserver(desMtx, wh_of_interest, onsets)
%
% Example:
% images = spm_get(Inf,'*img','Select images for timeseries plot');
% iv = InteractiveViewer('DisplayVol',images(1,:), 'AssociatedVols', images);
% desMtx = trialX;
% wh_of_interest = 1:144;
% attach(iv, tsFitPlotObserver(desMtx, wh_of_interest))
% % could also do: detach(iv, 1) then reattach
%

function obs = tsFitPlotObserver(desMtx, wh_of_interest, onsets)

    % Do things here that are not changed every time you click
    % Any variable I define here will be avail. in the inline, and thus in
    % the interactive viewer.
    
    [n, k] = size(desMtx);
    
    px = pinv(desMtx);
    
    desMtx_noint = desMtx;
    desMtx_noint(:, wh_of_interest) = [];
    
    px_noint = pinv(desMtx_noint);
    
    
    obs = IVObserver(@plot_);

    % This inline function is run whenever you click
    function plot_(iv, mmPos, voxPos)
        ts = get(iv, 'CurrentVolsTs');
        
        ndata = length(ts);
        if ndata < n
            disp('Warning! Timeseries has fewer data points than model!')
            desMtx = desMtx(1:ndata, :);
            px = px(:, 1:ndata);
            
            desMtx_noint = desMtx_noint(1:ndata, :);
            px_noint = px_noint(:, 1:ndata);
        end

        if ~iscol(ts), ts = ts'; end
        b = px * ts;
        
        % partial fit for columns of interest
        fit = desMtx(:, wh_of_interest) * b(wh_of_interest);
        
        % adjusted data, removing vars of no interest
        ts_adj = ts - desMtx_noint * px_noint * ts;
        
        hold off;
        plot(ts_adj, 'Color', 'k');
        
        hold on;
        plot(fit,'Color', 'r');
        
        
        title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f', mmPos(1), mmPos(2), mmPos(3)), 'FontSize', 14);
    end
end

