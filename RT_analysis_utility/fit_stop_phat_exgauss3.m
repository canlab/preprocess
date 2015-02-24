function [stop_phat,probstop_hat,sqerr_all] = fit_stop_phat_exgauss2(phat,ssoffset,rt,truestop)
% [stop_phat,probstop_hat] = fit_stop_phat_exgauss2(phat,ssoffset,rt,truestop)
%
% estimate the 2nd parameter (scale) of a gamma distribution for
% stop-signal processing times
% 
% minimize the distance between expected false alarms and observed false
% alarms, given 
%
            % we need an input phat w/one fixed param.  this is done above
%            
            % construct a function handle with fixed parameter
            % this will evaluate the pdf with one fixed argument
%            
% f is a function handle
% @(varparams) - these are output params that vary, i.e., are fitted
% f should return an error measure, e.g., sum of sq. errors, that you want
% to minimize
% create error measure: (est - true).^2
% objective function: thing to minimize: pstop from estimate_pstop, 
% which is estimate of % stops given parameters - true % stops
%f = @(phatstop) (estimate_pstop(rt,phat,phatstop,ssoffset) - truestop).^2;
%
% This version takes multiple rt distributions from different blocks as
% input.  do this once across all blocks!

f= @(phatstop) Calc_diff_exg_multiblk_orig(phatstop,phat,ssoffset,rt,truestop,'exgauss');

starttimes = [50:50:500];
for i = 1:length(starttimes),
    [sqerr(i)] = Calc_diff_exg_multiblk_orig(starttimes(i),phat,ssoffset,rt,truestop,'exgauss');
end
startvalue = starttimes(find(sqerr==min(sqerr)));
startvalue = startvalue(1);

% do 5 times at 3 different starting estimates; take median
starttimes = [startvalue-25,startvalue,startvalue+25];
for i = 1:length(starttimes), [stop_phat(i),fval(i)] = fminsearch(f,starttimes(i));, end

% this is the best-fitting values of phatstop
wh = find(fval==min(fval));       % stop latency (mean shift) in ms from stop sig onset
wh = wh(1);
stop_phat = stop_phat(wh);
sqerr_all = fval(wh);
% get the overall error and the fitted (est) probability of stopping for each
% block
[probstop_hat] = Calc_diff_exg_multiblk_orig(stop_phat,phat,ssoffset,rt,truestop,'exgauss');
            
return