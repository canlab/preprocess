function [gls,desc,prev] = rt_analyzer(y,varargin)
%
% input is y data (timeseries, designed primarily for RT data)
% followed by a series of cell arrays of strings containing
% string codes for 2-value factors
% 
% also handles autocorrelation, if flag 'ACF' is entered
% other flags are: 
% 'nback' to calculate regression and means
% based on previous trials
% 'prevcomp', followed by n, to look at y as a function of 
% percentage of trials in condition x in the last n trials
% 'nanexc' to exclude conditions with NaNs in means in desc
%
% 'accuracy' to get accuracy measures for each unique combo
% in desc.acc, and average only correct trials in desc.means
%
% gls contains beta weights regressing y against x values
%
% tor wager
%
%
% * assumes cell array inpus are cell arrays of strings 
%   with 2 unique non-NaN values and in column vectors
% 
% * observations are assumed to go forward in time down rows
%
% generate fake test data with:
% y = noisevector(300,[1 .4 .3 .2 .1],1);
% x(1:2:300) = {'on'}; x(2:2:300) = {'off'};
% x = x'; y = y';
% y = y+1; y(strcmp(x,'on')) = y(strcmp(x,'on'))+1;
% gls = rt_analyzer(y,x,'ACF');
% ...or...
%
% x(1:length(y)) = {'NaN'};, x(round(rand(150,1).*length(y)))={'on'};
% x(round(rand(150,1).*length(y)))={'off'}; x = x';
% gls = rt_analyzer(y,x,'ACF','nback',3);
%
% to get a null distribution:
% for i = 1:1000,
%   y = noisevector(300,[1 0],1)'; gls = rt_analyzer(y,x,'nback',3); b(i,:) = gls.betas;
% end
% for i = 1:length(gls.names),figure;hist(b(:,i),50),title(gls.names{i});,end
%
% or a distrib. where there should be an acf effect
% clear b
% for i = 1:1000,
%   y = noisevector(300,[1 .6 .3 0],1)'; gls = rt_analyzer(y,x,'nback',3,'ACF'); xc(i,:) = gls.xc;
% end
% for i = 1:length(gls.names),figure;hist(b(:,i),50),title(gls.names{i});,end
%
% More example code to get a distribution based on previous trial composition:
%
%[gls,desc,prev] = rt_analyzer(rts,Attsw,Obsw,'nback',1,'prevcomp',4); prev
%figure;plot(prev.xlevels,prev.trimmean(1:length(prev.count)./2),'rs-');
%hold on; plot(prev.xlevels,prev.trimmean(length(prev.count)./2+1:end),'bo-'); legend(prev.names);
%xlabel(['Previous trials ' prev.names{1} ' to ' prev.names{2}])
%
%figure;plot(prev(1).xlevels,prev(1).trimmean(1:length(prev(1).count)./2),'rs-');
%hold on; plot(prev(1).xlevels,prev(1).trimmean(length(prev(1).count)./2+1:end),'bo-'); legend(prev(1).names);
%figure;plot(prev(2).xlevels,prev(2).trimmean(1:length(prev(2).count)./2),'rs-');
%hold on; plot(prev(2).xlevels,prev(2).trimmean(length(prev(2).count)./2+1:end),'bo-'); legend(prev(2).names);



% ----------------------------------------
% * set up input arguments
% ----------------------------------------
ind = 1; nback = 0; doacf = 0; dointeract = 0; prevcomp = 0; nanexc = 0; doacc = 0;

for i = 1:length(varargin)
    if length(varargin{i}) > 1 & ~ischar(varargin{i}) % then an input
        inX{ind} = varargin{i}; 
        vname{ind} = inputname(i+1); 
        ind = ind+1;
        
    else
        switch varargin{i}
        case 'nback', nback = varargin{i+1};varargin{i+1} = 'done';
        case 'ACF', doacf = 1;
        case 'interact', dointeract = 1;
        case 'prevcomp', prevcomp = varargin{i+1};varargin{i+1} = 'done';
        case 'nanexc', nanexc = 1;
        case 'accuracy',acc = varargin{i+1}; doacc = 1; varargin{i+1} = 'done';
                
        end
    end
end




% ----------------------------------------
% * make columns
% ----------------------------------------
if isempty(inX), disp('No X data entered'),return, end

X = ones(size(inX{1}));                     % intercept; inX{1} is event name list (words)
u = {'intercept'}; 

tmp_names = unique(inX{1});
if length(tmp_names) == 1           % only 1 trial type
    x = X;                          % list of codes (1 -1), all 1s
    tmpu = tmp_names;               % 'contrast' -- here, just name
    descnms{1} = [tmp_names tmp_names tmp_names];            % names of descriptive stats (means), just name here
else
    [x,tmpu,descnms] = get_cols(inX{1},nback,dointeract);
end

for j = 1:length(tmpu), tmpu{j} = [vname{1} '_' tmpu{j}];,end
X = [X x]; u = [u tmpu];

for i = 2:length(inX), 
    [x,tmpu,tmpu2] = get_cols(inX{i},nback,dointeract);
    
    for j = 1:length(tmpu), tmpu{j} = [vname{i} '_' tmpu{j}];,end
    X = [X x]; u = [u tmpu]; descnms = [descnms tmpu2];
end

if size(X,1) ~= size(y,1), error('X and y lengths must be matched column vectors.'),end

whomit = find(isnan(y));
yo = y;
X(whomit,:) = [];
y(whomit) = [];

if doacc, acc(whomit) = [];,end
    
% ----------------------------------------
% * generalized least squares: 
%   acf model with pre-whitening
% ----------------------------------------
if doacf
    gls = genls(X,y);
else

% ----------------------------------------
% * fit model
% ----------------------------------------   
    gls.betas = (pinv(X) * y)';
    gls.X = X;
    gls.y = y;
end

gls.names = u;

% ----------------------------------------
% * Means of Conditions
%
% can check against this:
% tmp=inX{1}==1 & inX{2}==1 & inX{3}==1 & inX{4}==1 & inX{5}==0;
% tmp(whomit)=[];tmp=find(tmp);
% ----------------------------------------

% rows of c contain unique combinations of input indicator variables
% we want to find the mean RT for each unique combo of conditions

[c,wh] = unique(X,'rows'); %c = c(:,2:end); % remove intercept
if nanexc, c(any(c==0,2),:) = [];, end     % exclude rows with nan values

[desc.cnames] = build_names(c,u(2:end),descnms);

fstr = [repmat('%3.2f\t',1,size(c,1)) '\n'];
fstrn = [repmat('%3.0f\t',1,size(c,1)) '\n'];

for i = 1:size(c,1),
    % find all rows that match this row of c
    wh = ~any(repmat(c(i,:),size(X,1),1) - X,2);
    count(i) = sum(wh);
    
    % add accuracy, if entered
    if doacc, wh = wh & acc;,end
        
    wh = find(wh);    
    rt = y(wh); 
    
    if doacc, accuracy(i) = length(rt) ./ count(i);,end
    
    rt(isnan(rt)) = [];
    means(i) = mean(rt);
    % Now recursive trimming!!
    [tmp, ntrimmed(i)] = (trimts(rt,3,[])); %trimmean(rt,75);
    trimmeans(i) = mean(tmp);
    
    stdev(i) = std(tmp);
    sterr(i) = ste(tmp);
    
    
    % accuracy
    
end

desc.c = c;
desc.means = means;
desc.trimmeans = trimmeans;
desc.ntrimmed = ntrimmed;
desc.stdev = stdev;
desc.sterr = sterr;
desc.count = count;
if doacc, desc.acc = accuracy;,end
desc.str = sprintf(fstr,means);
desc.nstr = sprintf(fstrn,count);


% ----------------------------------------
% * Regression of y on previous trial 
%   composition n trials back
% ----------------------------------------
if prevcomp, 
    for i = 1:length(inX)
        prev(i) = get_prev(inX{i},prevcomp,yo,vname{i});,
    end
end

return






function [x,u,u2] = get_cols(X,nback,interact)

if iscell(X), 
    
    u = unique(X)'; u(strcmp(u,'NaN')) = [];
    x = zeros(size(X));
    if length(u) > 2, error('Only 2 unique non-NaN values in X allowed.'),end
    
    x(strcmp(X,u{1})) = 1;
    x(strcmp(X,u{2})) = -1;
    %x = x - mean(x);
    %x = x ./ 2;             % scale so unit effect is a - b
    
    % get names for descriptives
    %u2 = {u{2}}; if any(x==0), u2{2} = 'NaN';,end; u2 = [u2 {u{1}}];
    u2 = {[u(2) {'NaN'} u(1)]};
    
else        %if length(unique(X)) == 2   % it's a dummy code of some kind, leave alone
    u = unique(X)'; u(strcmp(u,'NaN')) = []; u(isnan(u)) = [];
    for i = 1:length(u), uu{i} = (u(i));, end, u = uu;   % u is cell array of strings
    x = zeros(size(X));
    if length(u) > 2, error('Only 2 unique non-NaN values in X allowed.'),end
    
    x(X == u{1}) = .5;
    x(X == u{2}) = -.5;
    
    for i = 1:length(u), u{i} = num2str(u{i});,end
    u2 = {[{(u{2})} {'NaN'} {(u{1})}]};
end

u = {[u{1} '-' u{2}]};

if nback > 0
    
    u2orig = u2;
    
    for i = 1:nback
        u2 = [u2 u2orig];
        
        clear xx
        xx(1:i,1) = NaN; %{'NaN'}; 
        xx = [xx; X];
        xx = xx(1:size(X,1));
        xx = get_cols(xx,0,0);
        x = [x xx]; u{end+1} = ['n-' num2str(i)];
    end
    
end

% interact would go here

return




function [names,nstr] = build_names(c,vname,descnms)

%rowstr = [repmat('%s\t',1,size(c,1)) '\n'];
names = [];
c = c(:,2:end); % remove intercept

for i = 1:size(c,2)         % rows in names
    
    for j = 1:size(c,1)     % cols in names
        
        if c(j,i) < 0, str = [vname{i} '_' descnms{i}{1}];
        elseif c(j,i) == 0, str = [vname{i} '_' descnms{i}{2}];
        elseif c(j,i) > 0, str = [vname{i} '_' descnms{i}{3}]; 
        else, error('this should never happen.'),
        end
        
        names = [names sprintf('%s\t',str)];
    end
    
    names = [names sprintf('\n')];
    
end

return



function out = get_prev(X,n,y,vname)
% X is cell array of cell arrays of string values denoting conditions

xprev(1:n,1) = NaN;

u = unique(X)'; u(strcmp(u,'NaN')) = [];
x = NaN .* zeros(size(X));
if length(u) > 2, error('Only 2 unique non-NaN values in X allowed.'),end
x(strcmp(X,u{1})) = -1; x(strcmp(X,u{2})) = 1;
    
for i = n+1:length(X)
    xprev(i) = mean(x(i-n:i-1));  
end
 
y(isnan(xprev)) = [];
x(isnan(xprev)) = [];
xprev(isnan(xprev)) = [];
xprev(isnan(y)) = [];
x(isnan(y)) = [];
y(isnan(y)) = [];

x = [ones(size(xprev)) x xprev xprev.^2 - mean(xprev.^2)];
out.bnames = {'Intercept' vname ['x_prev_' num2str(n)] 'x_prev^2'};
out.x = x;  

out.b = pinv(x) * y;
out.xlevels = unique(xprev);
out.names = u;
out.xnames = [' '];

for i = 1:length(out.xlevels)
    out.xmean(i) = mean(y(xprev == out.xlevels(i) & x(:,2) == -1));
    out.trimmean(i) = trimmean(y(xprev == out.xlevels(i) & x(:,2) == -1),75);
    out.count(i) = sum(xprev == out.xlevels(i) & x(:,2) == -1);
    out.xnames = sprintf('%s %s\t',out.xnames,[u{1} '_' num2str(out.xlevels(i))]);
    
end

for i = 1:length(out.xlevels)
    out.xmean(i+length(out.xlevels)) = mean(y(xprev == out.xlevels(i) & x(:,2) == 1));
    out.trimmean(i+length(out.xlevels)) = trimmean(y(xprev == out.xlevels(i) & x(:,2) == 1),75);
    out.count(i+length(out.xlevels)) = sum(xprev == out.xlevels(i) & x(:,2) == 1);
    out.xnames = sprintf('%s %s',out.xnames,[u{2} '_' num2str(out.xlevels(i))]);
end

out.str = ['levels run from ' u{1} ' to ' u{2} ', and within that from all prev. = ' u{1} ' (neg) to all ' u{2} ' (pos)'];

return






return
