choice = data{2}(:,1);
rew = data{2}(:,3);

cstick = zeros(size(choice, 1), 4); for i = 1:4, cstick(:, i) = choice == i; end
cstick = logical(cstick);

[t, decks] = size(cstick);

clear ev
% Expected value: Reward per selection
for i = 1:decks, ev(cstick(:,i), i) = rew(cstick(:,i)); end

decks_with_rewards = ev;

ev = cumsum(ev) ./ cumsum(cstick);  % expected value AFTER reward
tmp = NaN * zeros(1,4);  prior_ev = [tmp; ev(1:end - 1, :)];  % ev BEFORE each choice

prior_ev(isnan(prior_ev)) = 0;

tmp = prior_ev'; prior_ev_of_chosen = tmp(cstick');



% rank EV
rankev = rankdata(ev')';

rank_prior_ev = rankdata(prior_ev')'; %rank_prior_ev(any(isnan(prior_ev), 2), :) = 0;
tmp = rank_prior_ev'; rankprior_ev_of_chosen = tmp(cstick');


% switch vs. stay
stay = all(diff(cstick) == 0, 2);  % one indicates a stay on the subsequent trial
switch_deck = logical([~stay; 0]); % one indicates a switch on the next trial
stay = logical([stay; 0]);


p_switch = sum(switch_deck) ./ (length(switch_deck) - 1);
p_switch_after_neg = sum(switch_deck(rew < 0)) ./ length(switch_deck(rew < 0));
p_switch_after_pos = sum(switch_deck(rew > 0)) ./ length(switch_deck(rew > 0));

fprintf('\n%s\tp(switch)\tp(switch|rew)\tp(switch|pun)\t', filename);
fprintf('%3.2f\t%3.2f\t%3.2f\t', p_switch, p_switch_after_pos, p_switch_after_neg);

fprintf('Prior EV (chosen)\t%3.2f\tRank EV (chosen)\t%3.2f\n', nanmean(prior_ev_of_chosen), nanmean(rankprior_ev_of_chosen));



%%
f1 = create_figure('gambling plot'); figure(f1)

subplot(4, 1, 1)

hh = plot_onsets(rew .* (rew > 0), 'g', 1, 1);
hold on;
hh2 = plot_onsets(rew .* (rew < 0), 'r', 1, 1);

title('Rewards (green) and punishments (red)')
%%
subplot(4, 1, 2)

plot(prior_ev)
hold on;
ev2 = prior_ev;
ev2(~cstick) = NaN;
hold on; plot(ev2, '.-', 'LineWidth', 2)
title('Prior expected value (before feedback, chosen = dots), decks 1-4 = b, g, r, c')

%%

subplot(4, 1, 3)
plot(rank_prior_ev)
hold on;
rankev2 = rank_prior_ev;
rankev2(~cstick) = NaN;
hold on; plot(rankev2, '.-', 'LineWidth', 2)
title('Rank prior expected value (chosen = dots)')
set(gca, 'YLim', [0 5])


%% Transitions following rew/pun

rewstay = rew(stay);
rewswitch = rew(switch_deck);
nstay = length(rewstay);
nswitch = length(rewswitch);
xvals = ones(nstay, 1) + .1 * (rand(nstay, 1) - .5);

f1 = create_figure('choice plot'); figure(f1)

axh = gca; hold on; set(gcf, 'Color', 'w'); set(gca, 'FontSize', 16)
%axh = axes('Position', [.1 .07 .3 .2]);
hold on

plot([0 3], [0 0], 'Color', [.2 .2 .2]);

plot(xvals, rewstay, '.', 'Color', [.2 .2 .2], 'MarkerFaceColor', 'y')

xvals = 2 * ones(nswitch, 1) + .1 * (rand(nswitch, 1) - .5);
plot(xvals, rewswitch, '.', 'Color', [.2 .2 .2], 'MarkerFaceColor', 'y')

set(axh,'XLim', [0 3], 'XTick', [1 2], 'XTickLabel', {'Stay' 'Switch'});

ylabel('Reward preceding stay/switch');
title('Transitions following reward/pun');




