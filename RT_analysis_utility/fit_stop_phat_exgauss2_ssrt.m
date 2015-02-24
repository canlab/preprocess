function [stop_phat_ssrt,probstop_hat_ssrt,sqerr_ssrt,rsquare,avgdev,allrt] = fit_stop_phat_exgauss2_ssrt(phat,ssoffset,rt,truestop,options)
% [stop_phat,probstop_hat,sqerr,rsquare] = fit_stop_phat_exgauss2(phat,ssoffset,rt,truestop)
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

for i = 1:length(rt)
    rt{i}(isnan(rt{i})) = [];
end

f= @(ssrt) Calc_diff_exg_multiblk(ssrt,phat,ssoffset,rt,truestop,'ssrt');

% do 5 times at 3 different starting estimates; take median
starttimes = [50:50:500];
for i = 1:length(starttimes),
    [sqerr_ssrt(i)] = Calc_diff_exg_multiblk(starttimes(i),phat,ssoffset,rt,truestop,'ssrt');
end    
startvalue = starttimes(find(sqerr_ssrt==min(sqerr_ssrt)));
startvalue = startvalue(1);

starttimes = [startvalue-50 startvalue-25,startvalue,startvalue+25 startvalue+50];
for i = 1:length(starttimes), [stop_phat_ssrt(i),fval_ssrt(i)] = fminsearch(f,starttimes(i),options);, end

% this is the best-fitting values of phatstop
      % stop latency (mean shift) in ms from stop sig onset  
wh = find(fval_ssrt==min(fval_ssrt));
wh = wh(1);
stop_phat_ssrt = stop_phat_ssrt(wh); 
sqerr_ssrt = fval_ssrt(wh);
% get the overall error and the fitted (est) probability of stopping for each
% block

[sqerrterm, probstop_hat_ssrt] = Calc_diff_exg_multiblk(stop_phat_ssrt,phat,ssoffset,rt,truestop,'ssrt');

mse = sqerrterm ./ length(truestop);
% average deviation from true values for predicted values.
avgdev = sqrt(mse);

% prop. of variance across strategy blocks that's fitted (explained)
rsquare = (var(truestop) - var(truestop - probstop_hat_ssrt)) ./ var(truestop);

%sst = sum((truestop - mean(truestop)).^2);
%rsquare = 1 - (sqerr_ssrt./sst);

%plot_rt_hists_stopsig(rt,ssoffset,stop_phat_ssrt);


return