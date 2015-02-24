%%
% Establish batch loop
d = dir('*gamble_data.mat');

% all files need to have the same naming convention

%%
% Run Cell First to Set-Up Variables for Original Learning Portion of Task
% After running this cell, proceed to Individual Subject Analysis
disp('Original Learning Portion of Gambling Task');
appendname = 'gambling_original';
st = 1; en = 100;
doplot =0;
dosave = 0;
blocknames = {'Block1 Good' 'Block2' 'Block3' 'Block4' 'Block5'};
%%
% Run Cell to Set-Up Variables for REVERSAL  Portion of Task
% After running this cell, proceed to Individual Subject Analysis
disp('Reversal Learning Portion of Gambling Task');
appendname = 'gambling_reversal';
st = 101; en = 180;
doplot = 0;
dosave = 0;
blocknames = {'Block1 Good' 'Block2' 'Block3' 'Block4'};
%%
% Run Cell to Set-Up Variables for ALL TOGETHER  analysis
% After running this cell, proceed to Individual Subject Analysis
disp('Complete Gambling Task');
appendname = 'gambling_all';
st = 1; en = 180;
doplot = 1;
blocknames = {'Block1 Good' 'Block2' 'Block3' 'Block4' 'Block5' 'Block6' 'Block7' 'Block8' 'Block9'};
%%
dosave = 1;

%%
% Individual subject analysis
%
% Prints out a formatted table to be copied and pasted into excel. Also creates figures and saves to current working directory.
%
% Set st=1 and en= 100 and run script. Then set st=81 and en = 180 and run
% script again. Enter data side by side in table with first run labeled as
% "Original" and the 2nd one labeled as "Reversal"

group_stdrankchosen = zeros(en-st+1,length(d));
group_rank_v_chosen= zeros(en-st+1,length(d));
group_rank_ev_chosen = zeros(en-st+1,length(d));


nms = {'subid'};nms = [nms, blocknames];
nms = [nms,'Adv Small Picks' 'Adv Large' 'Dis Small' 'Dis Large' 'Avg Rank Chosen' 'Rank Hits' 'Modelfit' 'Loss Aversion' 'Decay Rate'];
nms = [nms,'Beta EV' 'Beta STD' 'Beta Freq' 'Rate Switchwin' 'Rate Switchloss'];
fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');

for i = 1:length(d),
    clear choice rew params best_params rank_error v pchoice negloglike rankv rankv_chosen rankev2 ev2 goodpicks winloss;
    load(d(i).name);
    choice = data{2}(st:en,1);
    rew = data{2}(st:en,3);
    fhan = @(params) gambling_exp_weighted_model(params, choice, rew);
    best_params = fminsearch(fhan, [1 .8]);
    [rank_error, v, pchoice, negloglike, rankv, rankv_chosen] = gambling_exp_weighted_model(best_params, choice, rew);
    avg_rankchosen = mean(rankv_chosen);
    chosetop = rankv_chosen==4;
    rank_hits = sum(chosetop);
    p_chance_rank = 1- binocdf(rank_hits,100,.25);

    good = strmatch('Advan',data{1});   %nominally good decks
    bad = strmatch('Disadvan',data{1}); %nominally bad decks

    goodpicks = zeros(en-st+1,1);
    for j = 1:length(goodpicks),
        goodpicks(j,1) = ismember(choice(j),good);
    end
    nblocks = length(goodpicks)./20;
    block_goodpicks = zeros(1,nblocks);
    for j = 0:nblocks-1,
        block_goodpicks(j+1) = nanmean(goodpicks(j*20+1:(j+1)*20));
    end
    
    switchwin = nan(length(goodpicks),1);
    switchloss = nan(length(goodpicks),1);

    actual_choices = zeros(length(goodpicks),4);

    wh1 = find(choice==1);
    wh2 = find(choice==2);
    wh3 = find(choice==3);
    wh4 = find(choice==4);

    actual_choices(wh1,1) = 1;
    actual_choices(wh2,2) = 1;
    actual_choices(wh3,3) = 1;
    actual_choices(wh4,4) = 1;

    wltmp = data{2}(st:en,3);

    winloss(wh1,1) = wltmp(wh1);
    winloss(wh2,2) = wltmp(wh2);
    winloss(wh3,3) = wltmp(wh3);
    winloss(wh4,4) = wltmp(wh4);

    for j = 1:(length(wltmp)-1),
        if winloss(j,choice(j)) > 0,
            if actual_choices(j+1,choice(j))==1,
                switchwin(j) = 0;
            else switchwin(j) = 1;
            end
        elseif winloss(j,choice(j)) < 0,
            if actual_choices(j+1,choice(j))==1,
                switchloss(j) = 0;
            else switchloss(j) = 1;
            end
        else disp('error');
        end
    end

    nswitchwin = length(switchwin)- sum(isnan(switchwin));
    nswitchloss = length(switchloss)- sum(isnan(switchloss));
    switchrate = (nswitchwin+nswitchloss)/(length(goodpicks)-1);
    rate_switchwin = nansum(switchwin)./nswitchwin;
    rate_switchloss = nansum(switchloss)./nswitchloss;
    
    

  
    
    
    ev = zeros(length(goodpicks),4);
    stdevrew = zeros(length(goodpicks),4);

    % ev is integrated across entire reward history
    % stdevrew is also integrated

    for j = 1:length(choice),
        window = sum(actual_choices(1:j,choice(j)));
        if j == 1,% first time step, initialize evs to 0
            ev(j,choice(j)) = winloss(j,choice(j));
            stdevrew(j,choice(j)) = 0;
        else
            ev(j,:) = ev(j-1,:);
            ev(j,choice(j)) = (sum(winloss(1:j,choice(j))))./sum(actual_choices(1:j,choice(j)));
            stdevrew(j,:) = stdevrew(j-1,:);
            tmp = winloss(1:j,choice(j)); tmp = tmp(tmp ~=0);
            if length(tmp) < j, newtmp = tmp(1:end); else newtmp = tmp(1:j); end
            stdevrew(j,choice(j)) = std(newtmp);
        end
    end

    prior_ev = [0 0 0 0; ev(1:end-1,:)];
    prior_v = [0 0 0 0; v(1:end-1,:)]; % model-based expected value

    rank_std = rankdata(stdevrew')';    % actual experienced stdev of rewards
    rank_priorev = rankdata(prior_ev')'; % actual experienced EV

    % get previous v (t - 1) for choices at time t
    prev_rankstd = cat(1, [2.5 2.5 2.5 2.5], rank_std(1:end-1, :));

%     rew = data{2}(:,3);
%     [rank_error, v, pchoice, negloglike, rankv, rankv_chosen] = gambling_exp_weighted_model(best_params, choice, rew);

    
    stdrankchosen = sum(prev_rankstd .* actual_choices,2);
    rank_priorev_chosen = sum(rank_priorev.*actual_choices,2);
    group_stdrankchosen(:,i) = stdrankchosen;  % not based on model
    group_rank_v_chosen(:,i) = rankv_chosen;   % based on model
    group_rank_priorev_chosen(:,i) = rank_priorev_chosen; % not based on model;

    
    % gets frequency of wins of each deck;
    [freqtmp,freqwins] = get_frequency_wins(data,st,en,actual_choices);
    
    x = [ev(st:en-1,1); ev(st:en-1,2);ev(st:en-1,3);ev(st:en-1,4)];
    x2 = [stdevrew(st:en-1,1); stdevrew(st:en-1,2);stdevrew(st:en-1,3);stdevrew(st:en-1,4)];
    x3 = [freqwins(st:en-1,1); freqwins(st:en-1,2); freqwins(st:en-1,3); freqwins(st:en-1,4)]; 
    y = [actual_choices(st+1:en,1); actual_choices(st+1:en,2); actual_choices(st+1:en,3);actual_choices(st+1:en,4)];
    [b,dev,stats] = glmfit([x,x2,x3],[y ones(size(y,1),1)],'binomial','logit'); 
    betas_EV = b(2);
    pvalues_logfit_EV = stats.p(2);
    betas_stdEV = b(3);
    pvalues_logfit_std = stats.p(3);
    betas_freqwin = b(4);
    pvalues_logfit_freqwin = stats.p(4);
    
    [Bev,dev,stats] = glmfit(x,[y ones(size(y,1),1)],'binomial','logit');
    [Bstd,dev,stats] = glmfit(x2,[y ones(size(y,1),1)],'binomial','logit');
    [Bfreq,dev,stats] = glmfit(x3,[y ones(size(y,1),1)],'binomial','logit');
    
    yfit = glmval(Bev,x,'logit');
    LL_ev = sum(log(yfit));
    yfit = glmval(Bstd,x2,'logit');
    LL_std = sum(log(yfit));
    yfit = glmval(Bfreq,x3,'logit');
    LL_freq = sum(log(yfit));
    
    reorder_decks(1) = strmatch('Advantageous large',data{1});
    reorder_decks(2) = strmatch('Advantageous small',data{1});
    reorder_decks(3) = strmatch('Disadvantageous large',data{1});
    reorder_decks(4) = strmatch('Disadvantageous small',data{1});

    if dosave
        fname = [d(i).name(1:4) '_' appendname];
        save(fname)
    end
    
%     if st == 101,
%         lastblock = nanmean(goodpicks(81:100)));
%     else lastblock = NaN;
%     end

    fprintf(1,'%s\t',d(i).name(1:4));
    fprintf(1,'%3.2f\t',block_goodpicks);
    fprintf(1,'%3.2f\t',sum(actual_choices(:,reorder_decks)));
    fprintf(1,'%3.2f\t',avg_rankchosen);
    fprintf(1,'%3.2f\t',rank_hits);
    fprintf(1,'%3.5f\t',p_chance_rank);
    fprintf(1,'%3.2f\t',best_params(1)); fprintf(1,'%3.2f\t',best_params(2));
    fprintf(1,'%3.5f\t',betas_EV);
    fprintf(1,'%3.5f\t',betas_stdEV);
    fprintf(1,'%3.5f\t',betas_freqwin);
    fprintf(1,'%3.5f\t',rate_switchwin);
    fprintf(1,'%3.5f\t',rate_switchloss);
    fprintf(1,'\n');

 
if doplot
    % plot figures

    %needs to be modified - plot of individual's choices as a function of EV of deck
    subplot(2, 1, 1)
    plot(prior_ev(:,reorder_decks))
    hold on;
    ev2 = prior_ev(:,reorder_decks);
    ev2(~actual_choices(:,reorder_decks)) = NaN;
    hold on; plot(ev2, '.-', 'LineWidth', 2)
    title('Prior nominal expected value (before feedback, not model based, chosen = dots), decks 1-4 = b, g, r, c')

    % needs to be modified
    % same plot but instead of raw EVs, ranks decks by EV
    subplot(2, 1, 2)
    plot(rankv(:,reorder_decks))
    hold on;
    rankev2 = rank_priorev(:,reorder_decks);
    rankev2(~actual_choices(:,reorder_decks)) = NaN;
    hold on; plot(rankev2, '.-', 'LineWidth', 2)
    title('Rank prior expected value (not model-based, chosen = dots)')
    set(gca, 'YLim', [0 5])

    if dosave
        name = [fname '_experienced_ev_choiceplots'];
        saveas(gcf,name, 'fig'); close(gcf);
    end

    figure; subplot(2, 1, 1)
    plot(prior_v(:,reorder_decks))
    hold on;
    ev2 = prior_v(:,reorder_decks);
    ev2(~actual_choices(:,reorder_decks)) = NaN;
    hold on; plot(ev2, '.-', 'LineWidth', 2)
    title('Prior nominal expected value (before feedback, model based, chosen = dots), decks 1-4 = b, g, r, c')

    subplot(2, 1, 2)
    plot(rankv(:,reorder_decks))
    hold on;
    rankev2 = rankv(:,reorder_decks);
    rankev2(~actual_choices(:,reorder_decks)) = NaN;
    hold on; plot(rankev2, '.-', 'LineWidth', 2);
    title('Rank prior expected value (model-based, chosen = dots)');
    set(gca, 'YLim', [0 5]);
    
    if dosave
        name = [fname '_modeled_ev_choiceplots'];
        saveas(gcf,name, 'fig'); close(gcf);
    end
end
end

save groupdata group*

%%
% Group Analysis

dat.groupnames = '1 = control, 2 = trier, 0 = not entered yet';
dat.model_rank_chosen = group_rankgoodev;
dat.rank_experienced_overall_std = group_stdrankchosen;

dat.mean_mrc = [mean(dat.model_rank_chosen(:, dat.group == 1), 2) mean(dat.model_rank_chosen(:, dat.group == 2), 2)];
figure; plot(dat.mean_mrc)
dat.mean_std = [mean(dat.rank_experienced_overall_std(:, dat.group == 1), 2) mean(dat.rank_experienced_overall_std(:, dat.group == 2), 2)];
figure; plot(dat.mean_std)
legend({'Control' 'Trier'})



cov = dat.group; cov(cov == 1) = -1; cov(cov == 2) = 1;
%h = timeseries_btwngroups_plot(dat.rank_experienced_overall_std', cov, [], 0);

dat_smooth = moving_average('gaussian', dat.rank_experienced_overall_std, 20);
h = timeseries_btwngroups_plot(dat_smooth', cov, [], 0);

plot_vertical_line(100);
legend(h, {'Control' 'Trier'});
scn_export_papersetup(500);
saveas(gcf,'std_of_chosen_by_group','png');

%%
% for selecting specific sections of timeseries for closer analysis, "wh"
% specifies trials to include
wh = [105:115];

ctrl_std = mean(dat.rank_experienced_overall_std(wh, dat.group == 1));
trier_std = mean(dat.rank_experienced_overall_std(wh, dat.group == 2));
ctrl_ev = mean(dat.model_rank_chosen(wh, dat.group == 1));
trier_ev = mean(dat.model_rank_chosen(wh, dat.group == 2));

create_figure('Bivariate plot');
plot(ctrl_ev, ctrl_std, 'bo','MarkerFaceColor', [.1 .1 .9]);
plot(trier_ev, trier_std, 'rs','MarkerFaceColor', [.7 .1 .1]);
xlabel('EV'); ylabel('STD')

