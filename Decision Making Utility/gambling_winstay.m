function [LMneg_winstay,pchosen] = gambling_winstay(data,st,en);
%[pchoice,LMneg,pchosen,actual_choices] = gambling_winstay(data)

actual_choices = zeros(100,4);      %
tmpchoice = data{2}(st:en,1);
wh1 = find(data{2}(st:en,1)==1);
wh2 = find(data{2}(st:en,1)==2);
wh3 = find(data{2}(st:en,1)==3);
wh4 = find(data{2}(st:en,1)==4);

actual_choices(wh1,1) = 1;
actual_choices(wh2,2) = 1;
actual_choices(wh3,3) = 1;
actual_choices(wh4,4) = 1;

deckprobs = (sum(actual_choices))./100;

winloss = data{2}(st:en,3);
choices = data{2}(st:en,1);
pchoice = zeros(100,4);

%%% I think what we need to do is estimate parameters of randonmness and
%%% sensitivity to outcomes, perhaps, and simulate the data under these
%%% parameters to see what simulations will give us the highest likelihood
%%% because individuals are not adopting entire win stay lose shift
%%% strategies and they're not all completely random

wh = find(actual_choices(1,:)>0);
pchoice(1,:) = 0; pchoice(1,wh) = 1;

for i = 2:100,
    
        if winloss(i-1) > 0,
            pchoice(i,:) = .0001;
            pchoice(i,choices(i-1)) = 1;
            % [1-myalpha base(2:end)./sum(base(2:end))*myalpha]
        elseif winloss(i-1) < 0,
            eligible_decks = find([1,2,3,4]~=choices(i-1));
            tmp = deckprobs(eligible_decks);
            tmp = tmp./sum(tmp);
            pchoice(i,eligible_decks) = tmp;
            pchoice(i,choices(i-1)) = .0001;
        else disp ('error')
        end
end

pchosen = sum(pchoice.*actual_choices,2);

LM = nansum(log(pchosen));
LMneg_winstay = -1*LM;

