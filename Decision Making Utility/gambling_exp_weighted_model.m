function [rank_error, v, pchoice, negloglike, rankv, rankv_chosen] = gambling_exp_weighted_model(params, choice, rew)
    %
    %
    %
    % [rank_error, v, pchoice, negloglike, rankv, rankv_chosen] = gambling_exp_weighted_model([1 .8], choice, rew);
    %
    % parameters:
    % [gam d], loss aversion (should be > 1 normally, 1 for no pref. weighting of losses) and history discounting (should be btwn 0 and 1)
    % history: 1 is no discounting.  so [1 1] should be the standard EV
    % model
    %
    % fhan = @(params) gambling_exp_weighted_model(params, choice, rew);
    % best_params = fminsearch(fhan, [1 .8])
    % [rank_error, v, pchoice, negloglike, rankv, rankv_chosen] = gambling_exp_weighted_model(best_params, choice, rew); rank_error
    % 
    % f1 = create_figure('Exponential p(choice)');
    % plot(pchoice)
    % pcchosen(~cstick) = NaN;
    % hold on; plot(pcchosen, '.-', 'LineWidth', 2)
    %
    % f1 = create_figure('Exponential p(choice)');
    % plot(rankv)
    % rankv_chosen(~cstick) = NaN;
    % hold on; plot(rankv_chosen, '.-', 'LineWidth', 2)

    % parameters
    gam = params(1);      % loss aversion; losses loom this much larger, try 1 or 1.2
    d = params(2);         % history discounting; value shrinks by this amount each trial, try .8; 1 = EV model


    cstick = zeros(size(choice, 1), 4); for i = 1:4, cstick(:, i) = choice == i; end
    cstick = logical(cstick);

    [T, decks] = size(cstick);

    % Expected value: Reward per selection
    for i = 1:decks, decks_with_rewards(cstick(:,i), i) = rew(cstick(:,i)); end

    clear r v

    for t = 1:T

        % current reward at time t, including loss aversion gam weighting
        % [r(t) - gam*p(t)]
        % -----------------------------------------------------
        r(t, :) = decks_with_rewards(t, :);
        r(t, r(t, :) < 0) = gam * r(t, r(t, :) < 0);

        % add weighted sum across prev rewards to get value, v
        % -----------------------------------------------------
        v(t, :) = r(t, :);
        denom = 1;

        for n = 1 : t-1

            v(t, :) = v(t, :) + d^n .* r(t - n, :);

            denom = denom + d^n;

        end

        % v(t) is value at time t, after reward feedback
        % should predict choice on NEXT trial
        % -----------------------------------------------------
        % divide by denominator, to track expected value rather than total
        % value
        v(t, :) = v(t, :) ./ denom;

    end

    % Summary measures

    % PCHOICE is NOT GOOD
    pchoice = (v.^2) ./ repmat(sum(v.^2, 2), 1, decks);

    pcchosen = pchoice;

    %%% idea: impose priors, based on overall values sampled from all
    %%% decks, and as information about a deck decreases, shrink towards
    %%% priors.  that way, info is imputed to un-chosen, or
    %%% non-recently-chosen, decks in a sensible way

    negloglike = -sum(log(pcchosen(cstick)));  % minimize this
    % issue: this is very sensitive to p-values near zero; thus, it's
    % dominated by the observations that fit poorly.

    rankv = rankdata(v')';
    
    % get previous v (t - 1) for choices at time t
    rankv = cat(1, [2.5 2.5 2.5 2.5], rankv(1:end-1, :));
    
    tmp = rankv'; rankv_chosen = tmp(cstick');

    rank_goodness = mean(rankv_chosen);
    rank_error = 4 - rank_goodness;  % minimize this.  will be sensitive to overall fit across observations,
    % but will lose something due to ranking


end




