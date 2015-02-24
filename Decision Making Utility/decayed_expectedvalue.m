function [pchoice,EVd,l,h,decay] = decayed_expectedvalue(a,data,st,en);

% Notes from Tor
% It'd be 0 for the current trial, plus sum(w(i)*r(i))
% for the previous trials where i = 1, 2, 3, etc... back

% this was about "objective" expected value, rather than subjective...
%     so r is the reward value on a given trial (just the number.)
%         

wh1 = find(data{2}(st:en,1)==1);
wh2 = find(data{2}(st:en,1)==2);
wh3 = find(data{2}(st:en,1)==3);
wh4 = find(data{2}(st:en,1)==4);


wltmp = data{2}(st:en,3); %wltmp = wltmp./1250;
winloss = zeros(length(wltmp),4);
winloss(wh1,1) = wltmp(wh1);
winloss(wh2,2) = wltmp(wh2);
winloss(wh3,3) = wltmp(wh3);
winloss(wh4,4) = wltmp(wh4);

EVd = zeros(length(wltmp),4);
choices = data{2}(st:en,1);
choicebytrial = {wh1, wh2, wh3, wh4};
decay = ones(length(wltmp),4);
l = zeros(length(wltmp),4);
h = zeros(length(wltmp),4);
objrew = zeros(length(wltmp),4); objrew(1,:) = winloss(1,:);
pchoice = zeros(length(choices),4);
firstchoices = [choicebytrial{1}(1),choicebytrial{2}(1),choicebytrial{3}(1),choicebytrial{4}(1)];

    
    
%sum(firstchoices) > 0% first time step, initialize evs to 0
    
for i=1:length(wltmp)
    if sum(i==firstchoices) > 0,
        deckfirst = find(i == [choicebytrial{1}(1),choicebytrial{2}(1),choicebytrial{3}(1),choicebytrial{4}(1)]);
        EVd(i,deckfirst) = winloss(i,deckfirst);
        if i > 1,
           for deck = 1:4,
            whc = find(choicebytrial{deck}<=i); %finds the index of the current value from the specific deck's list
            
            if length(whc) >1,
               whc = whc(end);
               l(i,deck) = choicebytrial{deck}(whc);
               h(i,deck) = i - l(i,deck);
               decay(i,deck) = (1+h(i,deck))^(-a);
               objrew(i,deck) = winloss(i,deck).*decay(i,deck);
               EVd(i,deck) = nansum(objrew(1:i,deck))./nansum(decay(1:i,deck));
            end
           end
        end
        
    else
        for deck = 1:4,
            whc = find(choicebytrial{deck}<=i); %finds the index of the current value from the specific deck's list
            
            if length(whc) >0,
               whc = whc(end);
               l(i,deck) = choicebytrial{deck}(whc);
               h(i,deck) = i - l(i,deck);
               decay(i,deck) = (1+h(i,deck))^(-a);
               objrew(i,deck) = winloss(i,deck).*decay(i,deck);
               EVd(i,deck) = sum(objrew(1:i,deck))./sum(decay(1:i,deck));
            end
        end

    end
end

mi = min(EVd(:)); ma = max(EVd(:));
decayed_ev_cent = (EVd-mi)./(ma-mi);

for i =2:length(wltmp),

 sumEVd = sum(decayed_ev_cent(i-1,:));
 pchoice(i,:) = decayed_ev_cent(i-1,:) ./sumEVd;
    
end



