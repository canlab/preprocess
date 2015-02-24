function [ev, ev_cent,stdevrew,freqwins,pchoice_ev,pchoice_std,pchoice_freq,LL_ev,LL_std,LL_freq] = getranks(data);
%[LMneg_EVrank,LMneg_varrank,
%betas_EV,pvalues_logfit_EV,betas_stdEV,pvalues_logfit_std] = getranks(data);

%%
actual_choices = zeros(100,4);      
tmpchoice = data{2}(1:100,1);
wh1 = find(data{2}(1:100,1)==1);
wh2 = find(data{2}(1:100,1)==2);
wh3 = find(data{2}(1:100,1)==3);
wh4 = find(data{2}(1:100,1)==4);

actual_choices(wh1,1) = 1;
actual_choices(wh2,2) = 1;
actual_choices(wh3,3) = 1;
actual_choices(wh4,4) = 1;

wltmp = data{2}(:,3); %wltmp = wltmp./1250;
winloss = zeros(length(wltmp),4);
winloss(wh1,1) = wltmp(wh1);
winloss(wh2,2) = wltmp(wh2);
winloss(wh3,3) = wltmp(wh3);
winloss(wh4,4) = wltmp(wh4);

ev = zeros(100,4);
stdevrew = zeros(100,4);
choices = data{2}(:,1);

for i = 1:100,
    window = sum(actual_choices(1:i,choices(i)));
    if i == 1,% first time step, initialize evs to 0
        ev(i,choices(i)) = winloss(i,choices(i));
        stdevrew(i,choices(i)) = 0;
    else
        ev(i,:) = ev(i-1,:);
        ev(i,choices(i)) = (sum(winloss(1:i,choices(i))))./sum(actual_choices(1:i,choices(i)));
        stdevrew(i,:) = stdevrew(i-1,:);
        tmp = winloss(1:i,choices(i)); tmp = tmp(tmp ~=0);
        if length(tmp) < i, newtmp = tmp(1:end); else newtmp = tmp(1:i); end
        stdevrew(i,choices(i)) = std(newtmp);
    end
end

%%

% myrank = zeros(100,4,nss);
% group_choices = zeros(100,4,nss);
% group_smoothed_choices = zeros(100,4,nss);
% bestrank_picks = zeros(100,nss);
% bestrank_smoothed = zeros(100,nss);
% indranks = zeros(100,4,nss);
% var_ranks = zeros(100,4,nss);
% rankvar = zeros(100,4,nss);
% reorder_decks = zeros(1,4);
% reordered_EVs = zeros(100,4,nss);
% reordered_stdev_EVs = zeros(100,4,nss);
% ste_ranks = zeros(100,4);
% ste_smchoices = zeros(100,4);
% EV_var_corr = zeros(nss,4);
% myrank_smoothed = zeros(100,4,nss);
% ste_myranks = zeros(100,4);
% rankvar_smoothed = zeros(100,4,nss);
% ste_rankvar = zeros(100,4);
% LMneg_EVrank = zeros(nss,1);
% LMneg_varrank = zeros(nss,1);
% LMneg_freqrank = zeros(nss,1);
% freq_ranks = zeros(100,4,nss);  %ranked(1-4) according to rate of wins (+ values)
% rankfreq = zeros(100,4,nss);
% ranksfreq_smoothed = zeros(100,4,nss);
% ste_rankfreq = zeros(100,4);

%%

% for ind = 1:length(sub_ind),
%     load (subdir(ind).name);
%     % make a new indexing based on the deck title to make sure everyone is
%     % in the same order
%     
%     %reorder expected values & variances according to the new deck order
%     reordered_EVs(:,:,ind) = group_EVs(:,reorder_decks,ind);
%     reordered_stdev_EVs(:,:,ind) = group_var_EV(:,reorder_decks,ind);
%     
% 
%     for i = 2:100,
%         indranks(i,:,ind) = rankdata(reordered_EVs(i-1,:,ind)'); % rank data according to the highest EV in new deck order
%         wh = find(actual_choices(i,reorder_decks) > 0); % find the actual choices based on this reordering
%         myrank(i,indranks(i,wh,ind),ind) = 1;   %myrank has columns in order of EV (first column is highest EV); 1 placed in the ranked column chosen by subjects;
%         var_ranks(i,:,ind) = rankdata(reordered_stdev_EVs(i-1,:,ind)');
%         rankvar(i,var_ranks(i,wh,ind),ind) = 1;
%     end
% 
% group_choices(:,:,ind) = actual_choices(:,reorder_decks);   % reorder the actual_choices to make the decks in the same order for everyone

% bestrank_picks(:,ind) = myrank(:,1,ind) + myrank(:,2,ind);
% bestrank_smoothed(:,ind) = moving_average('gaussian',bestrank_picks(:,ind),20);



% plot average ranks by expected value for all decks across trials
% (literally ranks 1-4, and plotting these values over time, but averaged
% across subjects)


%end % end loop through subject directory

%%
pchoice_ev = zeros(100,4);
pchoice_std = zeros(100,4);

mi = min(ev(:)); ma = max(ev(:));
ev_cent = (ev-mi)./(ma-mi);

 
for i=2:100,
    sumstr_ev = sum(ev_cent(i-1,:)); if sumstr_ev==0, sumstr_ev= 0.00001, end;
    pchoice_ev(i,:) = ev_cent(i-1,:) ./sumstr_ev;
    sumstr_std = sum(stdevrew(i-1,:)); if sumstr_std == 0, sumstr_std = 0.00001, end;
    pchoice_std(i,:) = stdevrew(i-1,:)./sumstr_std;
    %pchoice_std(i,:) = 1 - pchoice_std(i,:);
end
   
    pchosen_ev = sum(pchoice_ev .* actual_choices,2);
       wh = find(pchosen_ev==0);
       pchosen_ev(wh) = .0001;
       LM = nansum(log(pchosen_ev(20:100)));
       LMneg_EVrank = -1*LM;
    pchosen_std = sum(pchoice_std .* actual_choices,2);
       wh = find(pchosen_std==0);
       pchosen_std(wh) = .0001;
       LM = nansum(log(pchosen_std(20:100)));
       LMneg_varrank = -1*LM;

%% 
% figure; 
% t = 2:100;
% plot(t,meandeckrank(2:end,1),'b', t,meandeckrank(2:end,2),'g', t,meandeckrank(2:end,3),'r', t, meandeckrank(2:end,4),'c'); hold on;
% fill_around_line(meandeckrank(2:end,1),ste_ranks(2:end,1),'b'); fill_around_line(meandeckrank(2:end,2),ste_ranks(2:end,2),'g');
% fill_around_line(meandeckrank(2:end,3),ste_ranks(2:end,3),'r'); fill_around_line(meandeckrank(2:end,4),ste_ranks(2:end,4),'c');
% title('Actual EV Rank Across Time for Each Deck');
% xlabel('Trials 1-100');
% ylabel('Average Rank Based on Expected Value');
% legend('Advantageous small loss', 'Advantageous large loss', 'Disadvantageous small loss', 'Disadvantageous large loss');
% 
% 
% figure; 
% t = 2:100;
% plot(t,mean(myrank_smoothed(2:end,1,:),3),'b', t,mean(myrank_smoothed(2:end,2,:),3),'g', t, mean(myrank_smoothed(2:end,3,:),3),'r',t,mean(myrank_smoothed(2:end,4,:),3),'c'); hold on;
% fill_around_line(mean(myrank_smoothed(2:end,1,:),3),ste_myranks(2:end,1),'b'); fill_around_line(mean(myrank_smoothed(2:end,2,:),3),ste_myranks(2:end,2),'g');
% fill_around_line(mean(myrank_smoothed(2:end,3,:),3),ste_myranks(2:end,3),'r'); fill_around_line(mean(myrank_smoothed(2:end,4,:),3),ste_myranks(2:end,4),'c');
% title('Probability of Choosing From Best Ranked Decks');
% xlabel('Trials 1-100');
% ylabel('Probability of Choosing Deck');
% legend('Best Deck', '2nd Best Deck', '3rd Best Deck', 'Worst Deck');
% 
% figure; 
% t = 2:100;
% plot(t,mean(rankvar_smoothed(2:end,1,:),3),'b', t,mean(rankvar_smoothed(2:end,2,:),3),'g', t, mean(rankvar_smoothed(2:end,3,:),3),'r',t,mean(rankvar_smoothed(2:end,4,:),3),'c'); hold on;
% fill_around_line(mean(rankvar_smoothed(2:end,1,:),3),ste_rankvar(2:end,1),'b'); fill_around_line(mean(rankvar_smoothed(2:end,2,:),3),ste_rankvar(2:end,2),'g');
% fill_around_line(mean(rankvar_smoothed(2:end,3,:),3),ste_rankvar(2:end,3),'r'); fill_around_line(mean(rankvar_smoothed(2:end,4,:),3),ste_rankvar(2:end,4),'c');
% title('Probability of Choosing From Most Variable Deck');
% xlabel('Trials 1-100');
% ylabel('Probability of Choosing Deck');
% legend('Most Variable', '2nd in variance', '3rd in variance', 'Least variable');
% 
% 
% figure; 
% t = 2:100;
% plot(t,mean(group_smoothed_choices(2:end,1,:),3),'b', t,mean(group_smoothed_choices(2:end,2,:),3),'g', t, mean(group_smoothed_choices(2:end,3,:),3),'r',t,mean(group_smoothed_choices(2:end,4,:),3),'c'); hold on;
% fill_around_line(mean(group_smoothed_choices(2:end,1,:),3),ste_smchoices(2:end,1),'b'); fill_around_line(mean(group_smoothed_choices(2:end,2,:),3),ste_smchoices(2:end,2),'g');
% fill_around_line(mean(group_smoothed_choices(2:end,3,:),3),ste_smchoices(2:end,3),'r'); fill_around_line(mean(group_smoothed_choices(2:end,4,:),3),ste_smchoices(2:end,4),'c');
% title('Probability of Choosing From Each Deck Over Time');
% xlabel('Trials 1-100');
% ylabel('Probability of Choosing Deck');
% legend('Advantageous small loss', 'Advantageous large loss', 'Disadvantageous small loss', 'Disadvantageous large loss');

%%
[freqtmp,freqwins] = get_frequency_wins(data,1,100,actual_choices);
pchoice_freq = zeros(100,4);
    
%     reordered_choices= actual_choices(:,reorder_decks);
%     reordered_freqwins = freqwins(:,reorder_decks);
  
     for i = 2:100,
%         wh = find(reordered_choices(i,:) > 0);
%         freq_ranks(i,:) = rankdata(reordered_freqwins(i,:)');
%         deck_by_freq(i,freq_ranks(i,wh)) = 1;
            sumstr = sum(freqwins(i-1,:));
            pchoice_freq(i,:) = freqwins(i-1,:) ./sumstr;
     end


    
       pchosen = sum(pchoice_freq .* actual_choices,2);
       wh = find(pchosen==0);
       pchosen(wh) = .0001;
       LM = nansum(log(pchosen(20:100)));
       LMneg_freqrank = -1*LM;

       
%%
for i = 2:100,
        indranks(i,:) = rankdata(ev(i-1,:)'); % rank data according to the highest EV in new deck order
        wh = find(actual_choices(i,:) > 0); % find the actual choices based on this reordering
        myrank(i,indranks(i,wh)) = 1;   %myrank has columns in order of EV (first column is highest EV); 1 placed in the ranked column chosen by subjects;
end

rank1_last50 = (nansum(myrank(end-50+1:end,1)))./50;
rank2_last50 = (nansum(myrank(end-50+1:end,2)))./50;
rank1_last20 = (nansum(myrank(end-20+1:end,1)))./20;
rank2_last20 = (nansum(myrank(end-20+1:end,2)))./20;
%% logistic regression of EVs

    x = [ev(19:99,1); ev(19:99,2);ev(19:99,3);ev(19:99,4)];
    x2 = [stdevrew(19:99,1); stdevrew(19:99,2);stdevrew(19:99,3);stdevrew(19:99,4)];
    x3 = [freqwins(19:99,1); freqwins(19:99,2); freqwins(19:99,3); freqwins(19:99,4)]; 
    y = [actual_choices(20:100,1); actual_choices(20:100,2); actual_choices(20:100,3);actual_choices(20:100,4)];
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
%%
  fprintf(1,'\n'); fprintf(1,'\n');
   %nms = {'betaEV','pvalEV','betafreq', 'pvalfreq','betavar',%'pvalvar',...
        nms = {'LMnegEV', 'LMnegFreq', 'LMnegVar', 'betaEV','pvalEV','betafreq', 'pvalfreq','betavar', 'pvalvar',...
            'rank1_last50','rank2_last50','rank1_last20','rank2_last20'};
        fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');
        fprintf(1,'%3.4f\t',LMneg_EVrank); fprintf(1,'%3.4e\t',LMneg_freqrank); fprintf(1,'%3.4f\t',LMneg_varrank);
        fprintf(1,'%3.4f\t',betas_EV); fprintf(1,'%3.4f\t',pvalues_logfit_EV); fprintf(1,'%3.4f\t',betas_freqwin);
        fprintf(1,'%3.4f\t',pvalues_logfit_freqwin); fprintf(1,'%3.4f\t',betas_stdEV);fprintf(1,'%3.4f\t', pvalues_logfit_std);
        fprintf(1,'%3.4f\t',rank1_last50); fprintf(1,'%3.4f\t',rank2_last50);fprintf(1,'%3.4f\t',rank1_last20);fprintf(1,'%3.4f\t',rank2_last20);
        fprintf(1,'\n');
        fprintf(1,'\n');


%%
%win and loss sensitivity;
actual_choices = zeros(100,4);      

wh1 = find(data{2}(1:100,1)==1);
wh2 = find(data{2}(1:100,1)==2);
wh3 = find(data{2}(1:100,1)==3);
wh4 = find(data{2}(1:100,1)==4);

actual_choices(wh1,1) = 1;
actual_choices(wh2,2) = 1;
actual_choices(wh3,3) = 1;
actual_choices(wh4,4) = 1;

wltmp = data{2}(1:100,3); %wltmp = wltmp./1250;
rewards = cell(1,4);
switchdeck = cell(1,4);

winloss(wh1,1) = wltmp(wh1);   
winloss(wh2,2) = wltmp(wh2);
winloss(wh3,3) = wltmp(wh3);
winloss(wh4,4) = wltmp(wh4);

for deck = 1:4,
   for i =1:100,
       if winloss(i,deck) > 0 
           rewards{deck}(i,1) = winloss(i,deck);
           rewards{deck}(i,2) = 0;
       elseif winloss(i,deck) < 0,
           rewards{deck}(i,2) = winloss(i,deck);
           rewards{deck}(i,1) = 0;
       elseif winloss(i,deck) == 0,
           rewards{deck}(i,1) = 0;
           rewards{deck}(i,2) = 0;
       else disp('error');
       end
   end
end

switchwin = nan(100,1);
switchloss = nan(100,1);
choices = data{2}(:,1);

for i = 1:99,
  if winloss(i,choices(i)) > 0, 
    if actual_choices(i+1,choices(i))==1,
         switchwin(i) = 0;
    else switchwin(i) = 1;
    end
  elseif winloss(i,choices(i)) < 0,
      if actual_choices(i+1,choices(i))==1,
          switchloss(i) = 0;
      else switchloss(i) = 1;
      end
  else disp('error');
  end
end

nswitchwin = length(switchwin)- sum(isnan(switchwin));
nswitchloss = length(switchloss)- sum(isnan(switchloss));
switchrate = (nswitchwin+nswitchloss)/99;
rate_switchwin = nansum(switchwin)./nswitchwin;
rate_switchloss = nansum(switchloss)./nswitchloss;

 wins = [rewards{1}(19:99,1);  rewards{2}(19:99,1); rewards{3}(19:99,1); rewards{4}(19:99,1)]; 
 losses = [rewards{1}(19:99,2);  rewards{2}(19:99,2); rewards{3}(19:99,2); rewards{4}(19:99,2)]; 
 y = [actual_choices(20:100,1); actual_choices(20:100,2); actual_choices(20:100,3);actual_choices(20:100,4)];
 [b,dev,stats] = glmfit([wins,losses],[y ones(size(y,1),1)],'binomial','logit'); 
 beta_win = b(2);
 pvalues_win = stats.p(2);
 beta_loss = b(3);
 pvalues_loss = stats.p(3);
 
  fprintf(1,'\n'); fprintf(1,'\n');
   nms = {'betawin','pvalwin','betaloss','pvalloss','rate_switchwin','rate_switchloss' 'total switch rate'};
        fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');
        fprintf(1,'%3.4f\t',beta_win); fprintf(1,'%3.4f\t',pvalues_win); fprintf(1,'%3.4f\t',beta_loss); fprintf(1,'%3.4f\t',pvalues_loss); 
        fprintf(1,'%3.4f\t',rate_switchwin); fprintf(1,'%3.4f\t',rate_switchloss);fprintf(1,'%3.4f\t',switchrate;
        fprintf(1,'\n');
        fprintf(1,'\n');
