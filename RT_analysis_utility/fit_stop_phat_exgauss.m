function [stop_phatmu,stop_phatsig,probstop_hat] = fit_stop_phat_exgauss(rt,truestop,phatmix,phat,ssoffset)
% [stop_phat,probstop_hat] = fit_stop_phat_exgauss(rt,truestop,phatmix,phat,ssoffset)
%
% estimate the 1st parameter (mu) of an Ex-Gaussian distribution for
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
phatsig = 74;
f= @(phatstop) (Calc_diff_exgauss(phatstop,phatsig,0,ssoffset,rt)-truestop).^2;
        
% do 5 times at 3 different starting estimates; take median
starttimes = [100:50:400];
%vartimes = [0:50:300];
for i = 1:7, 
    [stop_phat,fval] = fminsearch(f,[starttimes(i)]);
    stop_phatmu(i) = stop_phat; %stop_phatsig(i) = stop_phat(2); 
    fval(i) = fval;
end
stop_phatmu = stop_phatmu(find(fval==min(fval)));       % stop latency (mean shift) in ms from stop sig onset
stop_phatmu = stop_phatmu(1); 
%stop_phatsig = stop_phatsig(find(fval==min(fval)));
stop_phatsig = phatsig;

%probstop_hat = estimate_pstop(rt,phat,stop_phat,ssoffset);
probstop_hat = Calc_diff_exgauss(stop_phatmu,stop_phatsig,0,ssoffset,rt);
            
return