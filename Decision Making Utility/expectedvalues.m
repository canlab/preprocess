function [group_EVs, group_var_EV] = expectedvalues(data,ind);

actual_choices = zeros(100,4);      
tmpchoice = data{2}(:,1);
wh1 = find(data{2}(:,1)==1);
wh2 = find(data{2}(:,1)==2);
wh3 = find(data{2}(:,1)==3);
wh4 = find(data{2}(:,1)==4);

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

% tmp1 = moving_average('gaussian',actual_choices(:,1),20); tmp1 = tmp1.*100;
% tmp2 = moving_average('gaussian',actual_choices(:,2),20); tmp2 = tmp2.*100;
% tmp3 = moving_average('gaussian',actual_choices(:,3),20); tmp3 = tmp3.*100;
% tmp4 = moving_average('gaussian',actual_choices(:,4),20); tmp4 = tmp4.*100;
% tmp5 = moving_average('gaussian',ev(:,1),20);
% tmp6 = moving_average('gaussian',ev(:,2),20);
% tmp7 = moving_average('gaussian',ev(:,3),20);
% tmp8 = moving_average('gaussian',ev(:,4),20);
% 
%         t=1:100;
%         figure; plot(t,tmp1,'b-',t,tmp5,'b-.',t,tmp2,'g-',t,tmp6,'g-.',...
%             t,tmp3,'r-',t,tmp7,'r-.',t,tmp4,'c-',t,tmp8,'c-.');
%         %title(files(1).name(1:3));
%         legend(data{1}{1}, 'EV', data{1}{2},'EV',...
%         data{1}{3},'EV',data{1}{4},'EV','Location','NW');

group_EVs(:,:,ind) = ev;
group_var_EV(:,:,ind) = stdevrew;
%group_freqwins(:,:,ind) = freqwins;
%group_choies(:,:,ind) = actual_choices;


clear LMneg a actual_choices ans bad c choices d data ev  ...
    f files fmode fval good goodpicks i k newtmp nms pchoice phats ...
    smgoodpicks start_a start_c start_w stdevrew strength t theta tmp ...
    tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmpchoice v w wh wh1 wh2 wh3 wh4 ...
    winloss wltmp x y y1 y2 z 

