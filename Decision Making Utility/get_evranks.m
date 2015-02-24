function [myrank,group_choices,rankvar,LMneg_EVrank,LMneg_varrank, EV_var_corr] = get_evranks(group_EVs,group_var_EV,subdir,sub_ind,nss);


myrank = zeros(100,4,nss);
group_choices = zeros(100,4,nss);
group_smoothed_choices = zeros(100,4,nss);
%bestrank_picks = zeros(100,nss);
%bestrank_smoothed = zeros(100,nss);
indranks = zeros(100,4,nss);
var_ranks = zeros(100,4,nss);
rankvar = zeros(100,4,nss);
reorder_decks = zeros(1,4);
reordered_EVs = zeros(100,4,nss);
reordered_stdev_EVs = zeros(100,4,nss);
ste_ranks = zeros(100,4);
ste_smchoices = zeros(100,4);
EV_var_corr = zeros(nss,4);
myrank_smoothed = zeros(100,4,nss);
ste_myranks = zeros(100,4);
rankvar_smoothed = zeros(100,4,nss);
ste_rankvar = zeros(100,4);
LMneg_EVrank = zeros(nss,1);
LMneg_varrank = zeros(nss,1);



for ind = 1:length(sub_ind),
    load (subdir(ind).name);
    % make a new indexing based on the deck title to make sure everyone is
    % in the same order
    
    %reorder expected values & variances according to the new deck order
    reordered_EVs(:,:,ind) = group_EVs(:,reorder_decks,ind);
    reordered_stdev_EVs(:,:,ind) = group_var_EV(:,reorder_decks,ind);
    

    for i = 2:100,
        indranks(i,:,ind) = rankdata(reordered_EVs(i-1,:,ind)'); % rank data according to the highest EV in new deck order
        wh = find(actual_choices(i,reorder_decks) > 0); % find the actual choices based on this reordering
        myrank(i,indranks(i,wh,ind),ind) = 1;   %myrank has columns in order of EV (first column is highest EV); 1 placed in the ranked column chosen by subjects;
        var_ranks(i,:,ind) = rankdata(reordered_stdev_EVs(i-1,:,ind)');
        rankvar(i,var_ranks(i,wh,ind),ind) = 1;
    end

group_choices(:,:,ind) = actual_choices(:,reorder_decks);   % reorder the actual_choices to make the decks in the same order for everyone

% bestrank_picks(:,ind) = myrank(:,1,ind) + myrank(:,2,ind);
% bestrank_smoothed(:,ind) = moving_average('gaussian',bestrank_picks(:,ind),20);



% plot average ranks by expected value for all decks across trials
% (literally ranks 1-4, and plotting these values over time, but averaged
% across subjects)


end % end loop through subject directory

for ind = 1:length(sub_ind),
    for col = 1:size(myrank,2),
        meandeckrank(:,col) = mean(indranks(:,col,:),3); % mean across all subjects
        group_smoothed_choices(:,col,ind) = moving_average('gaussian',group_choices(:,col,ind),20);
        EV_var_corr(ind,col) = corr(reordered_EVs(:,col,ind),reordered_stdev_EVs(:,col,ind)); %correlation between EV and variance per subject
        myrank_smoothed(:,col,ind) = moving_average('gaussian',myrank(:,col,ind),20); % smooth ranked picks of EV
        rankvar_smoothed(:,col,ind) = moving_average('gaussian',rankvar(:,col,ind),20);% smooth ranked picks of variance
    end
end


    for col = 1:size(myrank,2),
        for row = 2:size(myrank,1),
            for ind = 1:length(sub_ind),
                
            ste_ranks(row,col) = ste(indranks(row,col,:));
            ste_smchoices(row,col) = ste(group_smoothed_choices(row,col,:));
            ste_myranks(row,col) = ste(myrank_smoothed(row,col,:));
            ste_rankvar(row,col) = ste(rankvar_smoothed(row,col,:));
            end
        end
    end
    
for ind = 1:length(sub_ind),    
    pchosen = sum(myrank_smoothed(:,:,ind) .* group_choices(:,:,ind),2);    %smoothing doesn't make sense?
       wh = find(pchosen==0);
       pchosen(wh) = .0001;
       LM = nansum(log(pchosen));
       LMneg_EVrank(ind) = -1*LM;
    pchosen = sum(rankvar_smoothed(:,:,ind) .* group_choices(:,:,ind),2);
       wh = find(pchosen==0);
       pchosen(wh) = .0001;
       LM = nansum(log(pchosen));
       LMneg_varrank(ind) = -1*LM;
end


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




% logistic regression of EVs
betas_EV = zeros(nss,1);
betas_stdEV = zeros(nss,1);



pvalues_logfit_EV = zeros(nss,1);
pvalues_logfit_std = zeros(nss,1);


for ind = 1:length(sub_ind),
    x = [reordered_EVs(:,1,ind); reordered_EVs(:,2,ind);reordered_EVs(:,3,ind); reordered_EVs(:,4,ind)];
    x2 = [reordered_stdev_EVs(:,1,ind); reordered_stdev_EVs(:,2,ind);reordered_stdev_EVs(:,3,ind); reordered_stdev_EVs(:,4,ind)];
    y = [group_choices(:,1,ind); group_choices(:,2,ind); group_choices(:,3,ind); group_choices(:,4,ind)];
    [b,dev,stats] = glmfit([x,x2],[y ones(size(y,1),1)],'binomial','logit'); 
    betas_EV(ind) = b(2);
    pvalues_logfit_EV(ind) = stats.p(2);
    betas_stdEV(ind) = b(3);
    pvalues_logfit_std(ind) = stats.p(3);
    
end





