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
        ev(i,choices(i)) = (sum(winloss(i-window+1:i,choices(i))))./sum(actual_choices(i-window+1:i,choices(i)));
        stdevrew(i,:) = stdevrew(i-1,:);
        tmp = winloss(:,choices(i)); tmp = tmp(tmp ~=0);
        if length(tmp) < i, newtmp = tmp(1:end); else newtmp = tmp(1:i); end
        stdevrew(i,choices(i)) = std(newtmp);
    end
    
end


    for i = 2:100,
        indranks(i,:) = rankdata(ev(i-1,:)'); % rank data according to the highest EV in new deck order
        wh = find(actual_choices(i,:) > 0); % find the actual choices based on this reordering
        myrank(i,indranks(i,wh)) = 1;   %myrank has columns in order of EV (first column is highest EV); 1 placed in the ranked column chosen by subjects;
        %var_ranks(i,:,ind) = rankdata(reordered_stdev_EVs(i-1,:,ind)');
        %rankvar(i,var_ranks(i,wh,ind),ind) = 1;
    end

%group_choices(:,:,ind) = actual_choices(:,reorder_decks);   % reorder the actual_choices to make the decks in the same order for everyone

bestrank_picks1 = (nansum(myrank(:,1)))./99;
bestrank_picks2 = (nansum(myrank(:,2)))./99;

fprintf(1,'\n');
fprintf(1,'%3.4f\t',bestrank_picks1);
fprintf(1,'%3.4f\t',bestrank_picks2);
fprintf(1,'\n');
