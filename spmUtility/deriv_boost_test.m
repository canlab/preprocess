TR = 1;
secs = 500;
n = secs / TR; % images
n_events = 30; % number of events

res = 16;   % resolution, in samples per second


reps = 20;

true_intercept = 1;

amp = zeros(reps, 1);
intcpt = zeros(reps, 1);
betas = zeros(reps, 3);

boost = @(b) sign(b(1)) .* sqrt(sum(b .^ 2));

clear h t w

% NOTE:
% ONLY with orthonormal basis set, deriv boost amp seems to be stable
% across length & number of events

% GET HRF and NORMALIZE HRF the same way we normalized regressors, to get scaling factor
% ---------------------------------------------------------------------
bf = spm_get_bf(struct('name', 'hrf (with time and dispersion derivatives)', 'length', 30, 'dt', TR));
hrf = bf.bf;

% Mean-zero
hrf(:, 1:3) = scale(hrf(:, 1:3), 1);

% Orthogonalize Derivatives with respect to canonical, and 2nd wrt 1st
hrf(:, 1:3) = spm_orth(hrf(:, 1:3));

% Normalize for derivative boost
for ii = 1:3
    hrf(:, ii) = hrf(:, ii) ./ norm(hrf(:, ii));
end

scale_factor = max(hrf(:, 1)) .^ 2;
% NOT RIGHT, BUT I DON'T KNOW WHAT IS!

normed_hrf = hrf; % same stuff done to this as we're doing to the model mtx X
 
% ---------------------------------------------------------------------




create_figure('Example', 1, 3);
title('True Responses');

for r = 1:reps

    ybeta = rand(3, 1);         % [1; 0; 0]; %some arbitrary combo that produces an HRF-shaped response
    %ybeta(1) = .8 ; %ybeta(1) .* 2;  % Make shapes relatively more "canonical"

    % get onsets
    % ------------------------------------------
    ons = randperm(secs - 10);  % don't put events near the very end
    ons = ons(1:n_events);

    % get design matrix and basis functions
    % ------------------------------------------
    % norm option, if entered: create orthonbrmal basis set
    [X, delta, delta_hires, hrf] = onsets2fmridesign({ons'}, TR, secs, 'hrf (with time and dispersion derivatives)', 'norm');

    %[X, delta, delta_hires, hrf] = onsets2fmridesign({ons'}, TR, secs, 'hrf');

    
% %     % Mean-center
% %     X(:, 1:3) = scale(X(:, 1:3), 1);
% % 
% % 
% %     % Orthogonalize Derivatives with respect to canonical, and 2nd wrt 1st
% %     X(:, 1:3) = spm_orth(X(:, 1:3));
% % 
% %     % Normalize for derivative boost
% %     % Before model fitting
% % 
% %     for ii = 1:3
% %         X(:, ii) = X(:, ii) ./ norm(X(:, ii));
% %     end

    % Print check of norm, etc.
    if r == 1
        disp('Area:')
        areau = sum(X(:, 1:3) .^ 2)

        disp('Mean:')
        mean(X)

        disp('corr');
        disp(corrcoef(X(:, 1:3)));

    end



    % Create true response and "data" (outcome)
    % ------------------------------------------
    yresp = hrf * ybeta; % the "true" response

    yresp = yresp ./ max(yresp); % define max amplitude = 1

    subplot(1, 3, 1);
    plot(yresp, 'Color', rand(1,3)); drawnow

    y = true_intercept + conv(delta_hires, yresp); % the data (no noise, intercept = 1000, amp. = 1, LTI)
    y = y(1:res*n);

    % y is sampled at hi res (16 per TR), so downsample

    xi = (1:res:length(y))';       % downsample rate
    t = (1:length(y))';
    y = interp1(t, y, xi);



    % FIt and stuff

    b = pinv(X) * y;


    amp(r) = boost(b(1:3)) .* scale_factor;

    intcpt(r) = b(end);

    betas(r, :) = b(1:3)';


    subplot(1, 3, 2); hold off; plot(y); title('Example data time series'); hold on;
    plot(X * b, 'r-.')

    % IF WE RE-SCALE THE REGRESSORS, THIS IS GOING TO BE WRONG...FIX!
    % [h(r), t(r), w(r), auc(r)] = htw_from_fit(normed_hrf, b(1:3), TR, 'plot');
    [h(r), t(r), w(r), auc(r)] = htw_from_fit(hrf, b(1:3), TR ./ 16, 'plot', 'startval', 2.5, 'endval', 9);
pause(.5)

%     fit = normed_hrf * b(1:3) + b(end);
%     [h(r), t(r), w(r)] = fir2htw2(fit, 15 ./ TR, 0); % 15 sec / TR sampling resolution

end

% Plot
% ------------------------------------------


subplot(1, 3, 2); plot(y); title('Example data time series');
subplot(1, 3, 3); imagesc([scale(X(:, 1:end-1)) X(:, end)])
colormap gray; title('Example Design Matrix');

create_figure('Amplitude');
boxplot([betas amp h' intcpt]);
set(gca, 'XTick', 1:5, 'XTickLabel', {'B1' 'B2' 'B3' 'DB amplitude' 'Amp. from Fit' 'intercept'});

plot_horizontal_line(1)
text(1, 1.3, 'True amplitude', 'FontSize', 24);


% %             % Normalize for derivative boost
% %             areau = sum(hrf .^ 2);
% %             hrf = hrf ./ areau(ones(size(hrf, 1), 1), :);