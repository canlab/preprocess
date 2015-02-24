choices = zeros(100,1);
wltmp = zeros(100,1);
winloss = zeros(100,2);
v = zeros(100,1);
ev = zeros(100,4);
strength = ones(100,4);
pchoice = zeros(100,4);
wintmp = zeros(100,1);
losstmp = zeros(100,1);

w = .7;
a = .35;
c = .9;

deckindex = [3 3 3 3];

% deckvals1 = [50 -250 100 100 50 50 50 50 25 25]'; deckvals1 = getRandom(deckvals1);
% deckvals2 = [50 -50 -50 -50 -50 -50 150 100 100 100]'; deckvals2 = getRandom(deckvals2);
% deckvals3 = [100 -1250 150 150 150 100 100 100 100 50]'; deckvals3 = getRandom(deckvals3);
% deckvals4 = [100 -150 -200 -250 -300 -350 250 250 200 200]'; deckvals4 = getRandom(deckvals4);
deckwins1 = [50 50 50 50 50 50 50 50 50 50]'; %deckwins1 = getRandom(deckwins1);
deckloss1 = [0 0 -50 0 -50 0 -50 0 -50 -50]'; %deckloss1 = getRandom(deckloss1);
deckwins2 = [50 50 50 50 50 50 50 50 50 50]'; %deckwins2 = getRandom(deckwins2);
deckloss2 = [0 0 0 0 0 0 0 0 0 -250]'; %deckloss1 = getRandom(deckloss1);
deckwins3 = [100 100 100 100 100 100 100 100 100 100]'; %deckwins1 = getRandom(deckwins1);
deckloss3 = [0 0 -150 0 -300 0 -200 0 -250 -350]'; %deckloss1 = getRandom(deckloss1);
deckwins4 = [100 100 100 100 100 100 100 100 100 100]'; %deckwins2 = getRandom(deckwins2);
deckloss4 = [0 0 0 0 0 0 0 0 0 -1250]'; %deckloss1 = getRandom(deckloss1);


%cardarray = [deckvals1, deckvals2, deckvals3, deckvals4];
cardwins = {deckwins1 deckwins2 deckwins3 deckwins4};
cardloss = {deckloss1 deckloss2 deckloss3 deckloss4};

choices(1:8) = [1;2;3;4;1;2;3;4];
% wltmp(1:8) = [cardarray(1,choices(1)); cardarray(1,choices(2)); cardarray(1,choices(3)); cardarray(1,choices(4));...
%     cardarray(2,choices(1)); cardarray(2,choices(2)); cardarray(2,choices(3)); cardarray(2,choices(4))];
wintmp(1:8) = [cardwins{1}(1);cardwins{2}(1);cardwins{3}(1);cardwins{4}(1);cardwins{1}(2); cardwins{2}(2);cardwins{3}(2);cardwins{4}(2)];
losstmp(1:8) = [cardloss{1}(1);cardloss{2}(1);cardloss{3}(1);cardloss{4}(1);cardloss{1}(2); cardloss{2}(2);cardloss{3}(2);cardloss{4}(2)];

t = 1:length(winloss);          % trial numbers
theta = (t/10).^c; 

for i = 1:8,
%     if wltmp(i) > 0,
%         winloss(i,1) = wltmp(i);
%         winloss(i,2) = 0;
%     elseif wltmp(i) <0,
%         winloss(i,1) = 0;
%         winloss(i,2) = wltmp(i);
%     end
    winloss(i,:) = [wintmp(i),losstmp(i)];

    v = w .* winloss(:,1) + (1 - w) .* winloss(:,2);
    
    if i == 1,% first time step, initialize evs to 0
        ev(i,choices(i)) = a.*v(i);
    else
        ev(i,:) = ev(i-1,:);
    ev(i,choices(i)) = a.*v(i) + (1-a).*ev(i-1, choices(i));
    end
end

 for j = 8:100,
     
%      if wltmp(j) > 0,
%         winloss(j,1) = wltmp(j);
%         winloss(j,2) = 0;
%     elseif wltmp(j) <0,
%         winloss(j,1) = 0;
%         winloss(j,2) = wltmp(j);
%     end
    winloss(j,:) = [wintmp(j),losstmp(j)];
    v = w .* winloss(:,1) + (1 - w) .* winloss(:,2); 

    
    ev(j,:) = ev(j-1,:);
    ev(j,choices(j)) = a.*v(j) + (1-a).*ev(j-1, choices(j));
    

    if j == 8,
        strength(j,1) = exp(theta(j).*ev(j,1)); strength(j,2) = exp(theta(j).*ev(j,2)); strength(j,3) = exp(theta(j).*ev(j,3)); strength(j,4) = exp(theta(j).*ev(j,4));
    elseif j > 8,
        strength(j,:) = strength(j-1,:);
        strength(j,choices(j)) = exp(theta(j).*ev(j,choices(j)));
        %strength(j,choices(j)) = min([7.5667e+304,strength(j,choices(j))]);
        strength(j,choices(j)) = max([0.00000001,strength(j,choices(j))]);
        %strength(j,choices(j)) = max([0.00000001,strength(j,choices(j))]);
    end

    if length(find(choices(1:j)==choices(j))) >=40
        strength(j,choices(j)) = 0;
    end

    
    
    sumstr = sum(strength(j,:));
    pchoice(j+1,:) = strength(j,:) ./sumstr;

    wh = find(pchoice(j+1,:)==max(pchoice(j+1,:)));
    wh = wh(1);
    %randval = rand;
%     if randval < .99,
%         choices(j+1) = wh;
%     else 
%         tmp = [1;2;3;4];
%         tmp = getRandom(tmp);
%         choices(j+1) = tmp(1);
%     end
    choices(j+1) = wh;
    deckindex(choices(j+1)) = deckindex(choices(j+1)) + 1;
    if deckindex(choices(j+1)) > 10,
        deckindex(choices(j+1)) = 1;
    end
    %wltmp(j+1) = cardarray(deckindex(choices(j+1)),choices(j+1));
    wintmp(j+1) = cardwins{choices(j+1)}(deckindex(choices(j+1)));
    losstmp(j+1) = cardloss{choices(j+1)}(deckindex(choices(j+1)));
    
 end
 
 wh1 = find(choices==1);
 wh2 = find(choices==2);
 wh3 = find(choices==3);
 wh4 = find(choices==4);
 
actual_choices(wh1,1) = 1;
actual_choices(wh2,2) = 1;
actual_choices(wh3,3) = 1;
actual_choices(wh4,4) = 1;

choices = choices(1:100);
actual_choices = actual_choices(1:100,:);
%wltmp = wltmp(1:100);
wintmp = wintmp(1:100);
losstmp = losstmp(1:100);

data = cell(1,4);
data{2}(:,1) = choices;
pchoice = pchoice(1:100,:);
pchosen_sim = sum(pchoice .* actual_choices,2);
data{2}(:,3) = wltmp;
%data{4} = pchosen_sim;
data{2}(:,3) = wintmp;
data{2}(:,4) = losstmp;
