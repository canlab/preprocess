function [pchoice,v,ev,strength,theta] = predict_choices(phats,data,st,en)

%[pchoice,v,ev,strength,theta] = predict_choices([w,a,c],data)
% s is a vector of one Sj value for each deck, e.g., [.8 .2 .3 .15]
% returns pchoice, a vector of prob. of choices for each of the 4 decks


% 
%dat = {deck_names, trial_data, total_score}
%deck names = data{1}{:}
%trial_data = {deckchoice, rt, trialscore}
%deck choice = data{2}(:,1);
%rt = data{2}(:,2);
%trialscore = data{2}(:,3);
phats(1) = max([0 phats(1)]);   % phat(1) = w, and we want to constrain it to be between 0 and 1
phats(1) = min([1 phats(1)]);
phats(2) = max([0 phats(2)]);
phats(2) = min([1 phats(2)]);
% phats(3) = max([0 phats(3)]);

w = phats(1);
a = phats(2);
c = phats(3);


% choices = [1 0 0 0; 0 1 0 0; etc.]
% winloss = [10 0; 0 -5; 20 0; etc.]
%wltmp = data{2}(st:en,3); %wltmp = wltmp./1250;
wintmp = data{2}(st:en,3); losstmp = data{2}(st:en,4);
winloss = ones(length(wintmp),2);
% for i=1:length(wltmp),
%     if wltmp(i) > 0,
%         winloss(i,1) = wltmp(i);
%         winloss(i,2) = 0;
%     elseif wltmp(i) <0,
%         winloss(i,1) = 0;
%         winloss(i,2) = wltmp(i);
%     else fprintf(1,'error','s');
%     end
% end
winloss = [wintmp losstmp];
%choiceindex = [1 4 3 2 4 3 2 3 4]

t = 1:length(winloss);          % trial numbers

theta = (t/10).^c;           % sensitivity for each trial
%theta = repmat(c,100);
%theta = c.^t;
% note: different in diff papers ***

% ---------------- Get Valence -----------------
%
% ----------------------------------------------
choices = data{2}(st:en,1);
% for i = 1:length(t),
%     wh = choices(i);
%     if i ==1,
%         v(i,:) = 0;
%         v(i,wh) = w.*winloss(i,2) + (1-w).*winloss(i,1);
%     else
%         v(i,:) = v(i-1,:);
%         v(i,wh) = w.*winloss(i,2) + (1-w).*winloss(i,1);
%     end
% end

%v = get_valence(w,winloss);     % this returns valence for ALL trials
v = w .* winloss(:,1) + (1 - w) .* winloss(:,2);   %totally wrong



% ---------------- Get expected Values for each deck ---------------
%
% ------------------------------------------------------------------

ev = zeros(length(choices),4);      % zero expected value for each deck

    % get expected value for each deck;
           % which deck they chose
strength = ones(length(choices),4);

for i = 1:length(t)
    if i == 1,% first time step, initialize evs to 0
        ev(i,choices(i)) = a.*v(i);
    else
        ev(i,:) = ev(i-1,:);
        ev(i,choices(i)) = a.*v(i) + (1-a).*ev(i-1, choices(i)); % expectancy for each deck
        %ev(i,choices(i)) = ev(i-1,choices(i)) +
        %%a.*(v(i)-ev(i-1,choices(i)));    % equation according to 2006 paper
    end
end


% ------------------ Get strength of Decks -------------------------
%
% ------------------------------------------------------------------


for i = 1:length(t)
    if i == 1,
        %%expectancy valuence model approach with exponential term in
        %%strength equation
        %K = 500; % constant to avoid negatives;
        strength(i,choices(i)) = exp(theta(i).*ev(i,choices(i)));
        %strength(i,choices(i)) = theta(i).*ev(i,choices(i)); %modification
    else
        strength(i,:) = strength(i-1,:);
        strength(i,choices(i)) = exp(theta(i).*ev(i,choices(i)));
        %strength(i,choices(i)) = theta(i).*ev(i,choices(i));
        strength(i,choices(i)) = min([7.5667e+304,strength(i,choices(i))]);
        strength(i,choices(i)) = max([0.00000001,strength(i,choices(i))]);
    end
    
    if length(find(choices(1:i)==choices(i))) >= 40
           strength(i,choices(i)) = 0;
    end
%     wh = choices(i);
%     strength(i,:) = strength(i-1,:);
%     strength(i,wh) = exp(theta(i).*ev(i,wh));           % strength for each deck 
 
  
end

pchoice = zeros(length(choices),4);

warning off
for i = 2:length(t),
    sumstr = sum(strength(i-1,:));
    pchoice(i,:) = strength(i-1,:) ./sumstr;
end
warning on

%pchoice = pchoice(1:100,:);
    



%return