function o2 = spm_results_display_canlab(t_image_name, p_thresh, varargin)
% Display results from a standard SPM analysis
%
% o2 = spm_results_display_canlab(t_image_name, p_thresh, [meth])
%
% meth = 'orthviews', 'montage', or 'all' [default]
%
% o2 = spm_results_display_canlab('montage', 'spmT_0003.img', .001);

load SPM

% create statistic_image object.
% p-values automatically calculated because we have entered 'type' 't' and dfe.
t = statistic_image('image_names', t_image_name, 'type', 't', 'dfe', SPM.xX.erdf);

% Threshold and display
t = threshold(t, p_thresh, 'unc');

if length(varargin) > 0
    meth = varargin{1};
end

switch meth
    
    case {'all', 'orthviews'}
        
        orthviews(t)
        
        % To change the colormap:
        cm = spm_orthviews_change_colormap([0 0 1], [1 1 0], [0 .5 1], [0 .7 .5], [1 0 .5], [1 .5 0]);
        
    case {'all', 'montage'}
        
        r = region(t);
        
        o2 = canlab_results_fmridisplay(r, 'noblobs', 'nooutline');
        
        o2 = addblobs(o2, r, 'splitcolor', {[0 0 1] [0 .5 1] [1 0 .5] [1 1 0]});
         
    otherwise
        error('Unknown display method')
        
end
