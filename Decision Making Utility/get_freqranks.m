%function [LMneg_freqrank,betas_freqwin,pvalues_logfit_freqwin] = get_freqranks(data,st,en);


freq_ranks = zeros(100,4);  %ranked(1-4) according to rate of wins (+ values)
deck_by_freq = zeros(100,4);
ste_rankfreq = zeros(100,4);
actual_choices = zeros(100,4);      

tmpchoice = data{2}(st:en,1);
wh1 = find(data{2}(st:en,1)==1);
wh2 = find(data{2}(st:en,1)==2);
wh3 = find(data{2}(st:en,1)==3);
wh4 = find(data{2}(st:en,1)==4);

actual_choices(wh1,1) = 1;
actual_choices(wh2,2) = 1;
actual_choices(wh3,3) = 1;
actual_choices(wh4,4) = 1;


    reorder_decks(1) = strmatch('Advantageous small',data{1});
    reorder_decks(2) = strmatch('Advantageous large',data{1});
    reorder_decks(3) = strmatch('Disadvantageous small',data{1});
    reorder_decks(4) = strmatch('Disadvantageous large',data{1});
    
    [freqtmp,freqwins] = get_frequency_wins(data,st,en,actual_choices);
    
%     reordered_choices= actual_choices(:,reorder_decks);
%     reordered_freqwins = freqwins(:,reorder_decks);
%   
     for i = 2:100,
%         wh = find(reordered_choices(i,:) > 0);
%         freq_ranks(i,:) = rankdata(reordered_freqwins(i,:)');
%         deck_by_freq(i,freq_ranks(i,wh)) = 1;
            sumstr = sum(freqwins(i-1,:));
            pchoice(i,:) = freqwins(i-1,:) ./sumstr;
     end


    
       pchosen = sum(pchoice .* actual_choices,2);
       wh = find(pchosen==0);
       pchosen(wh) = .0001;
       LM = nansum(log(pchosen(10:100)));
       LMneg_freqrank = -1*LM;


%logistic regression of current win history on choices



    x = [freqwins(11:99,1); freqwins(11:99,2); freqwins(11:99,3); freqwins(11:99,4)]; 
    y = [actual_choices(12:100,1); actual_choices(12:100,2); actual_choices(12:100,3); actual_choices(12:100,4)];
    [b,dev,stats] = glmfit(x,[y ones(size(y,1),1)],'binomial','logit'); 
    betas_freqwin = b(2);
    pvalues_logfit_freqwin = stats.p(2);


   fprintf(1,'\n'); fprintf(1,'\n');
        nms = {'LMnegEV', 'LMnegFreq', 'LMnegVar', 'betaEV','pvalEV','betafreq', 'pvalfreq','betavar', 'pvalvar'};
        fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');
        fprintf(1,'%3.4f\t',LMneg_EVrank); fprintf(1,'%3.4e\t',LMneg_freqrank); fprintf(1,'%3.4f\t',LMneg_varrank);
        fprintf(1,'%3.4f\t',betas_EV); fprintf(1,'%3.4f\t',pvalues_logfit_EV); fprintf(1,'%3.4f\t',betas_freqwin);
        fprintf(1,'%3.4f\t',pvalues_logfit_freqwin); fprintf(1,'%3.4f\t',betas_stdEV);fprintf(1,'%3.4f\t', pvalues_logfit_std);
        fprintf(1,'\n');
        fprintf(1,'\n');

% figure; 
% t = 2:100;
% plot(t,mean(rankfreq_smoothed(2:end,1,:),3),'b', t,mean(rankfreq_smoothed(2:end,2,:),3),'g', t, mean(rankfreq_smoothed(2:end,3,:),3),'r',t,mean(rankfreq_smoothed(2:end,4,:),3),'c'); hold on;
% fill_around_line(mean(rankfreq_smoothed(2:end,1,:),3),ste_rankfreq(2:end,1),'b'); fill_around_line(mean(rankfreq_smoothed(2:end,2,:),3),ste_rankfreq(2:end,2),'g');
% fill_around_line(mean(rankfreq_smoothed(2:end,3,:),3),ste_rankfreq(2:end,3),'r'); fill_around_line(mean(rankfreq_smoothed(2:end,4,:),3),ste_rankfreq(2:end,4),'c');
% title('Probability of Choosing From Deck with Most Frequent Wins');
% xlabel('Trials 1-100');
% ylabel('Probability of Choosing Deck');
% legend('Advantageous small loss', 'Advantageous large loss','Disadvantageous small loss', 'Disadvantageous large loss');