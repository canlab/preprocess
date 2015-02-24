function [ncrit,pow,obspow,str,ncrit50] = power_calc(SNR,alph,varargin)
% function [ncrit,pow,obspow,str,ncrit50] = power_calc(SNR (cohen's d) t-value or z,alpha or tcrit,[n observed],[d, t, z, or r flag; assumes d],[plot color 'k--'])
% 
% calcuates number of subjects required to achieve significance
% for a one-sample t-test
%
% Tor Wager
%
% Both correlation and d inputs tested and work, Jan 30, 06
%
% INPUTS:
%
% SNR signal to noise ratio, est / stdev of estimate
%   or t value (but enter n for t-stat input as the last argument)
% 
% tcrit  alpha value, one-tailed (e.g., .025) or critical t given
% correction for multiple comparisons; one-tailed
%
% n number of subjects, if t statistic
% tcrit critical t value for significance
%
% OUTPUT: critical number of subjects given this effect size.
% -- pow is power for n = 2:1000
% --if optional n is entered, also gives observed power given n
% Error if obs n > 1000 !!
%
% Tor Wager
%
% Examples:
% [ncrit,pow,obspow] = power_calc(.4,.025,21); ncrit,obspow
% [ncrit,pow,obspow] = power_calc(.4,.025,21,'d','k--'); ncrit,obspow
% [ncrit,pow,obspow] = power_calc(.4,.025,21,'r','k--');    % correlation
% of .4

% SNR = d / se; CV = se / d;  (zcrit*CV)^2

%% t = b / (std ./ sqrt(n)) = b * sqrt(n) / std
% t / sqrt(n) = b / std = SNR

% t + tinv(.2,??df) > 0

%eff = .5; sigma = 1; n = 30; df = n-1; 
%t = eff./(sigma ./ sqrt(n))
%tcrit = tinv(.975,df)
%beta = tcdf(tcrit - t,df), pow = 1 - beta

SNR = abs(SNR);

if length(varargin) > 0,
    nobs = varargin{1};
else
    nobs = 20; % pick an arbitrary value; affects obspow, but not pow and ncrit
end

% CONVERT everything to d effect size measure 
if length(varargin) > 1,
    switch varargin{2}
        
    case 't',
        % SNR is a t-value
        SNR = SNR ./ sqrt(nobs);    % cohen's d, the SNR, in original sample sd units
    case 'z'
        %z-score
        z = SNR;    
        p = normcdf(z);
        t = tinv(p,nobs - 1);
        SNR = t ./ sqrt(nobs);
    
    case 'r'
        % correlation effect size
        r = SNR;
        %r = SNR ./ ((SNR.^2 + 4).^.5); % d to r, from cohen (1988)
        % simulate r:
        % for i = 1:5000, x = mvnrnd(ones(19,2),[1 .4472; .4472 1]); [a,pp]
        % = corrcoef(x); p(i) =pp(1,2); r(i) = a(1,2);,end
        
        SNR = (2*r) / ((1 - r.^2).^.5);
            
    case 'd'
        % cohen's d    
        % do nothing
        
    otherwise 
        error('Flag must be ''t'' or ''d'' or ''z''');
    end
end

if alph < 1 % then it's an alpha value (e.g., .025) -- nothing to do
    %tcrit = tinv(1-alph,nobs - 1);
else
    % this is a critical t-value, convert to alph
    alph = 1 - tcdf(alph,nobs - 1);
end

% notes:
% http://userwww.sfsu.edu/~efc/classes/biol710/power/power%20analysis%20web
% .htm#4%20power%20of%20test
% says: d = t * sqrt(na+nb / nanb), or t*sqrt(2/n) for two-sample equal n.
% checked by hand.  t-values are sqrt(2) smaller for equal var 2-sample t.

% generate power for all n
n = 2:500;
df = n - 1;
t = SNR .* sqrt(n);


if exist('r') == 1, 
    % do this for correlation
    t = (r.*sqrt(n-2))/(sqrt(1-r.^2));
end

tcrit = tinv(1-alph,df);
beta = tcdf(tcrit - t,df);
pow = 1 - beta;

ncrit = n(find(beta <= .2));
if isempty(ncrit),
    ncrit = 'Inf';
else
    ncrit = ncrit(1);
end

ncrit50 = n(find(beta <= .5));
if isempty(ncrit50),
    ncrit50 = 'Inf';
else
    ncrit50 = ncrit50(1);
end


if length(varargin) > 0,
    obspow = pow(find(nobs == n));
else
    obspow = [];
end

str = sprintf('One sample: d = %3.2f, N needed 80%% = %3.0f , N needed 50%% = %3.0f\n',SNR,ncrit, ncrit50);
disp(str)


% two-group, equal N, pooled sigma
t = t ./ sqrt(2);
df = n - 2;
beta = tcdf(tcrit - t,df);
ncrit_2g = n(find(beta <= .2));
ncrit_2g50 = n(find(beta <= .5));

if isempty(ncrit_2g),
    ncrit_2g = 'Inf';
else
    ncrit_2g = ncrit_2g(1);
end
if isempty(ncrit_2g50),
    ncrit_2g50 = 'Inf';
else
    ncrit_2g50 = ncrit_2g50(1);
end

str = sprintf('2-group: d = %3.2f, N per group 80%% = %3.0f , N per group 50%% = %3.0f\n',SNR,ncrit_2g, ncrit_2g50);
disp(str);    



doplot = []; % doplot is color, if not 0
if length(varargin) > 2, doplot = varargin{3};,end
if doplot, 
    hold on; set(gca,'FontSize',16)
    plot(n,pow,doplot,'LineWidth',3);,
    %plot([0 1000],[.8 .8],'k-');
    %h = text(ncrit+5,.8,str); set(h,'FontSize',16);
    
    % add to legend
    [LEGH,OBJH,OUTH,OUTM] = legend;
    str = {[input('Enter name for this line: ','s') ' ' str]};
    if ~isempty(OUTM), OUTM(end+1) = str;, else, OUTM = str;,end
    legend(OUTM);

end


return




    % get power, given t-value and alpha
    %tcrit = tinv(.975,n - 1);
    %beta = tcdf(tcrit - SNR,n - 1);, 
    %pow = 1 - beta; % power for this test
    
    % now get critical n for a power of .8
    %beta = .2;  % beta = 1-pow = 1-.8
    
    %SNR = SNR ./ sqrt(n);  
    
    
