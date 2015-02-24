function [freqtmp,freqwins] = get_frequency_wins(data,st,en,actual_choices);

wltmp = data{2}(st:en,3); %wltmp = wltmp./1250;
winloss = zeros(length(wltmp),4);
wh1 = find(data{2}(st:en,1)==1);
wh2 = find(data{2}(st:en,1)==2);
wh3 = find(data{2}(st:en,1)==3);
wh4 = find(data{2}(st:en,1)==4);
winloss(wh1,1) = wltmp(wh1);
winloss(wh2,2) = wltmp(wh2);
winloss(wh3,3) = wltmp(wh3);
winloss(wh4,4) = wltmp(wh4);

freqtmp = zeros(100,4);
freqwins = zeros(100,4);
choices = data{2}(:,1);

for i = 1:100,
    window = sum(actual_choices(1:i,choices(i)));
    if i == 1,% first time step, initialize evs to 0
        if winloss(i,choices(i)) > 0,
            freqtmp(i,choices(i)) = 1; 
        else freqtmp(i,choices(i)) = 0;
        end
        
        wh = find(choices(1:i) == choices(i)); whfreq = freqtmp(wh,choices(i));
        freqwins(i,choices(i)) = (sum(whfreq(end-window+1:end)))./length(whfreq);
    else
        freqwins(i,:) = freqwins(i-1,:);
        if winloss(i,choices(i)) > 0,
            freqtmp(i,choices(i)) = 1; 
        else freqtmp(i,choices(i)) = 0;
        end
        wh = find(choices(1:i) == choices(i));whfreq = freqtmp(wh,choices(i));
        freqwins(i,choices(i)) = (sum(whfreq(end-window+1:end)))./length(whfreq);
        %freqwins(i,choices(i)) = (sum(freqtmp(i-window+1:i,choices(i))))./sum(actual_choices(i-window+1:i,choices(i)));
    end
    
end

%group_freqwins(:,:,ind) = freqwins;