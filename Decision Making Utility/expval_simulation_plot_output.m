% notes: july 9
% strength is unreasonable for all choices of c
% pchoice is perfect, even for sens. (c) = 0
% model sometimes chooses invalid things?  try c = .5.  should always pick
% highest EV?
% theta and theta * trial could just be free params

figure; subplot(4, 1, 1)

% % plot(choices); set(gca, 'YLim', [0 5])
% % title('Choices')

%%
subplot(4, 1, 1)

hh = plot_onsets(winloss(:,1), 'g', 1, 1);
hold on;
hh2 = plot_onsets(winloss(:,2), 'r', 1, 1);
plot(v)
title('Valence (blue), wins (green), and losses (red)')


subplot(4, 1, 2)
plot(ev);  
ev2 = ev;
cstick = zeros(size(choices, 1), 4); for i = 1:4, cstick(:, i) = choices == i; end
ev2(~cstick) = NaN;
hold on; plot(ev2, 'LineWidth', 2)

hold on; plot(ev2, '.-', 'LineWidth', 2)
title('Expected value, EV')

%%
subplot(4, 1, 3)

plot(strength);  
strength2 = strength;
strength2(~cstick) = NaN;

hold on; plot(strength2, '.-', 'LineWidth', 2)
title('Strength')

%%
subplot(4, 1, 4)

plot(pchoice);  
pchoice2 = pchoice;
pchoice2(~cstick) = NaN;

hold on; plot(pchoice2, '.-', 'LineWidth', 2)
title('P Choice')

