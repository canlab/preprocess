function [stop_phat,probstop_hat] = fit_stop_phat(rt,truestop,phat,ssoffset)
% [stop_phat,probstop_hat] = fit_stop_phat(rt,truestop,phat,ssoffset)
%
% estimate the 2nd parameter (scale) of a gamma distribution for
% stop-signal processing times
% 
% minimize the distance between expected false alarms and observed false
% alarms, given 

            % we need an input phat w/one fixed param.  this is done above
            
            % construct a function handle with fixed parameter
            % this will evaluate the pdf with one fixed argument
            
% f is a function handle
% @(varparams) - these are output params that vary, i.e., are fitted
% f should return an error measure, e.g., sum of sq. errors, that you want
% to minimize
% create error measure: (est - true).^2
% objective function: thing to minimize: pstop from estimate_pstop, 
% which is estimate of % stops given parameters - true % stops
%f = @(phatstop) (estimate_pstop(rt,phat,phatstop,ssoffset) - truestop).^2;
f= @(phatstop) (Calc_diff_exgauss(phat(1),phat(2),phat(3),phatmix(i,1),phatmix(i,2),phatmix(i,3))-truestop).^2;
            
% do 5 times at 3 different starting estimates; take median
for i = 1:5, stop_phat(i) = fminsearch(f,15 * i);, end
stop_phat = median(stop_phat(i));
            
probstop_hat = estimate_pstop(rt,phat,stop_phat,ssoffset);

pstop(ssr,r) = Calc_diff_exguass(phat(1),phat(2),phat(3),varystop(ssr),phat(2),phat(3));
            
return