function [stop_phatmu,stop_phatsig,probstop_hat,sqerr_all,rsquare,avgdev] = fit_stop_phat_exgauss2(phat,ssoffset,rt,truestop,options)
% [stop_phat,probstop_hat] = fit_stop_phat_exgauss2(phat,ssoffset,rt,truestop)
%
% estimate the 1st parameter(mu) of ex-Gaussian distribution for
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

f= @(phatstop) Calc_diff_exg_multiblk(phatstop(1),phatstop(2),ssoffset,rt,truestop,'exgauss');

fprintf(1,'Sampling search space...');
starttimes = [50:50:500];
vartimes = [25:25:250];
for i = 1:length(starttimes),
    for k=1:length(vartimes),
    [sqerr(i,k),tmp] = Calc_diff_exg_multiblk(starttimes(i),vartimes(k),ssoffset,rt,truestop,'exgauss');
    end
end
[r,c] = matrix_min(sqerr);
startvaluemu = starttimes(r);
startsigma = vartimes(c);

fprintf(1,'Nonlinear least squares');
% do 5 times at 3 different starting estimates; take median
starttimes = [startvaluemu-50,startvaluemu-25,startvaluemu,startvaluemu+25,startvaluemu+50];
vartimes = [startsigma-50,startsigma-25,startsigma,startsigma+25,startsigma+50];
warning off
for i = 1:length(starttimes),
    fprintf(1,'\b\b\b%03d',i);
    for k = 1:length(vartimes),
    [stop_phat,fval(i,k)] = fminsearch(f,[starttimes(i), vartimes(k)],options); 
    stop_phatmu(i,k) = stop_phat(1); stop_phatsig(i,k) = stop_phat(2); 
    end
end
warning on
    
% this is the best-fitting values of phatstop
[r,c] = matrix_min(fval);
stop_phatmu = stop_phatmu(r,c);
stop_phatsig = stop_phatsig(r,c);
sqerr_all = fval(r,c);
% get the overall error and the fitted (est) probability of stopping for each
% block
[sqerrterm, probstop_hat] = Calc_diff_exg_multiblk(stop_phatmu,stop_phatsig,ssoffset,rt,truestop,'exgauss');

mse = sqerr_all ./ length(truestop);
% average deviation from true values for predicted values.
avgdev = sqrt(mse);

% prop. of variance across strategy blocks that's fitted (explained)
rsquare = (var(truestop) - var(truestop - probstop_hat)) ./ var(truestop);

%sst = sum((truestop - mean(truestop)).^2);
%rsquare = 1 - (sqerr_ssrt./sst);

%plot_rt_hists_stopsig(mixedplots,rt,ssoffset,stop_phat,'EI predicted stops:');
            
return

function [r,c] = matrix_min(a)

minimum = min(a(:));
minimum = find(a == min(a(:)));

if length(minimum) > 1, warning('Multiple exact matches for min.');
    minimum = minimum(1);, end

z = zeros(size(a)); z(minimum) = 1;
[r,c] = find(z);

return