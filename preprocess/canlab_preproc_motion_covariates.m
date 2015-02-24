function [mvmtcat, mvmt_by_session, motion_cov_set] = canlab_preproc_motion_covariates(mvmt_param_files, images_per_session, varargin)
% [mvmtcat, mvmt_by_session, motion_cov_set] = canlab_preproc_motion_covariates(mvmt_param_files, images_per_session, varargin)
%
% Create motion covariates and plot them if requested.
% used in canlab_preproc_2012.m
%
% motion_cov_set: 'Movement params plus ^2, diffs, ^2 diffs (zscored), sep sessions';
% For use as nuisance covariates in linear models.
%
% To plot:
% [...] = canlab_preproc_motion_covariates(PREPROC.mvmt_param_files, images_per_session, 'plot')
%
% To plot and save:
% [...] = canlab_preproc_motion_covariates(PREPROC.mvmt_param_files, images_per_session, 'savedir', pwd)
%         * Enter dir name instead of pwd to choose save location.
%
% In canlab_preproc:
% [~, ~, PREPROC.Nuisance.motion_cov_set] = canlab_preproc_motion_covariates(PREPROC.mvmt_param_files, PREPROC.images_per_session, 'savedir', PREPROC.qcdir);
%
% Tor Wager, Aug 2012
%
% Examples:
% mvmtfiles = filenames('remi1495/Functional/rp_a*txt')
% [~, ~, motion_cov_set] = canlab_preproc_motion_covariates(mvmtfiles, [660 330 660 330]);

doplot = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'savedir'
                savedir = varargin{i+1};
                varargin{i+1} = [];
                doplot = 1;
                
            case 'plot', doplot = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Load
% -------------------------------------------------------------------------

for i = 1:length(mvmt_param_files)
    
    if exist(mvmt_param_files{i}, 'file')
        mvmt{i} = load(mvmt_param_files{i});
        
    else
        warning('Warning!!! Movement parameter file missing.', mvmt_param_files{i});
    end
end

mvmtcat = vertcat(mvmt{:});

mvmt_by_session = blkdiag(mvmt{:});

% create cov set from higher-order transformations
% -------------------------------------------------------------------------
for i = 1:length(mvmt)
    sessm{i} = zscore(mvmt{i});
    
    smpd = [0 0 0 0 0 0; diff(sessm{i})];
    
    sessmplus{i} = [sessm{i} sessm{i}.^2 smpd smpd.^2];
end

nuisance.sessmvmtplus = blkdiag(sessmplus{:});

motion_cov_set = blkdiag(sessmplus{:});

% Plot
% -------------------------------------------------------------------------

if doplot
    
    create_figure('Movement covariates', 2, 1)
    
    subplot(2, 1, 1)
    m = mvmtcat(:, 1:3);
    draw_session_boxes(images_per_session, m(:))
    h = plot(m);
    axis tight;
    xlabel('Time'); ylabel('mm');
    legend(h, {'x' 'y' 'z'});
    
    subplot(2, 1, 2)
    m = mvmtcat(:, 4:6);
    draw_session_boxes(images_per_session, m(:))
    h = plot(m);
    axis tight;
    xlabel('Time'); ylabel('Radians');
    legend(h, {'roll' 'pitch' 'yaw'});
    
    sz = get(0, 'ScreenSize');
    set(gcf, 'Position', [50 50 sz(3) ./ 2 sz(4) ./ 3]);
    
    if exist('savedir', 'var')
        
        if ~exist(savedir, 'dir'), mkdir(savedir), end
        scn_export_papersetup(600);
        saveas(gcf, fullfile(savedir, 'movement_analysis.png'));
    end
    
end

end % main function


% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------


function draw_session_boxes(images_per_session, y)

if diff(size(images_per_session)) < 0, images_per_session = images_per_session'; end

nsess = length(images_per_session);
colors = {'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k'};  %{'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
while nsess > length(colors), colors = [colors colors]; end

yval = min(y(:)) - .05 * min(y(:));
yheight = (max(y(:)) + .05 * max(y(:))) - yval;

c = cumsum(images_per_session);
st = [0 c(1:end-1)];

for jj = 1:2:nsess
    h1 = drawbox(st(jj), images_per_session(jj), yval, yheight, colors{jj}(1));
    set(h1, 'FaceAlpha', .10, 'EdgeColor', 'none');
end

end

