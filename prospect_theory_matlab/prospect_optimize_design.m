%% setup

xvals = [25 50 75 100 150 200 400 800];
yvals = [-800 -400 -200 -150 -100:25:100 150 200 400 800];
xvals = yvals;
pvals = [.01:.02:1];
pop = combvec(xvals,yvals,pvals)';
wh = find(pop(:,1)==pop(:,2));
pop(wh,:) = [];
x = pop(:,1);
y = pop(:,2);
p = pop(:,3);

truev = .7;
truea = .5;
trueb = .5;
truelo = 1.5;
truep = [truea trueb truev truelo];

iter = 100;
%u = (x.^truev) .* wfcn(p,truea,trueb) + ((y .^ truev) .* (1-wfcn(p,truea,trueb))) + 5 .* randn(length(u),1);
%u is true u in head space.  put through inverse value function to go back
%to monetary value scale.

ntrials = 50;
ntotalgamb = size(pop,1);
gindx = 1:ntotalgamb;

%% Try out intermediate commands (just to show you...)

% select a random subset
wh = randsample(gindx,ntrials,'true');

% get fitness
fitness = prospect_organism(wh,pop,truep,iter);

% alternatively, using function handles:
genfun = @() randsample(gindx,ntrials,'true')';
wh = genfun();

objhan = @(wh) prospect_organism(ceil(wh),pop,truep,iter);
objhan(wh);


%% Run GA

iter = 100; %iterations per organism to get st. errs
numgen = 100;
gensize = 100;

objhan = @(wh) prospect_organism(ceil(wh),pop,truep,iter);
genfun = @() randsample(gindx,ntrials,'true')';

% starting input; only size matters, not content
wh = randsample(gindx,ntrials,'true');

[best_params,fit,beff,in] = tor_ga(gensize,numgen,wh,objhan,genfun);

xyp = pop(best_params{1},:);

save ga_finished best_params fit beff in pop ntrials