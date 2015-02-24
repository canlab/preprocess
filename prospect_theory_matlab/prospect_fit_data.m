function [best_params,fit_u,fval] = prospect_fit_data(value1,y,p,rep_price,doplot)
% best_params = prospect_fit_data(x,y,p,rep_price);
% 
% x = [25 50 75 100 150 200 400 800 50 75 100 150 150 200 200]; x = [x x];
% y = [0  0  0  0   0   0   0   0   25 50 50  50  100 100 150]; y = [y -y];
% p = rand(size(x));
% u = x .* p + (y .* (1-p));

 % x = best outcome
    % y = alternative outcome
    % p = prob. of getting x
    
    % a = scale, weighting function
    % b = exp., weighting function
    % m = exp., value function
    % lo = loss aversion
    
    
% x = data(:,1);
% y = data(:,2);
% p = data(:,3);
% u = data(:,4);

if nargin < 5, doplot = 0; end

startparams = [.5 .5 .7 2];

% params = cell(4,1);
% params{1}(1,1) = 0.001; params{1}(:,:,2) = 4;
% params{2}(1,1) = 0.001; params{2}(:,:,2) = 4;
% params{3}(1,1) = 0.001; params{3}(:,:,2) = 4;
% params{4}(1,1) = 0.001; params{4}(:,:,2) = 4;


% create a handle to the link function, to pass into error function
%fhan = @(params,x,y,p) prospect_utility(params,x,y,p);
%fhan = @(params,value1,y,p) uinv(prospect_utility_mixedgambles(params,value1,y,p),params(3));
%fhan = @(params,value1,y,p) prospect_utility_mixedgambles(params,value1,y,p);
fhan = @(params,value1,y,p) temp_prospectfit(params,value1,y,p);;

% objective: SSE of the vinv of utility, to transform back to orig. scale
% units like reported price is
% equation 4
%objfun = @(params) sum((rep_price - fhan(params,value1,y,p).^(1/params(3))).^2);
objfun = @(params) sum((rep_price - fhan(params,value1,y,p)).^2);

% if we could assess true underlying utility equivalent, it would be:
%objfun = @(params) sum((rep_price - fhan(params,x,y,p)).^2);



%options = optimset('Display','off','MaxIter',100000,'LevenbergMarquardt','on');

    %PROBLEM = struct('objective',objfun,'x0',startparams,'options',options,'solver','fminsearch');
warning off
    [best_params fval] = fminsearch(objfun,startparams);
    %[best_params,fit,beff,in,isconverged] = tor_ga(81,50,params,objfun);
    warning on
    
    
    if nargout > 1
        fit_u = fhan(best_params,value1,y,p);

    end
    
    
    if doplot
        xvals = [-100:100]; len = length(xvals);
        u = prospect_utility(best_params,xvals,zeros(1,len),.5 * ones(1,len));
        %tor_fig(1,2);
        subplot(1,2,1);
        plot(xvals,u,'LineWidth',2);
        plot_vertical_line(0);
        title('Utility Function (at p = .5)')
        xlabel('Nominal value');
        ylabel('Est. Subjective Utility')

        p = [.001:.001:.999];
        wp = wfcn(p,best_params(1),best_params(2));
        subplot(1,2,2);
        plot(p,wp,'LineWidth',2)
        title('Probability Weighting Function')
        xlabel('Probability');
        ylabel('Subjective Probability')
        drawnow
    end
    
end




function wp = wfcn(p,a,b)
    % function wp = wfcn(p,a,b)
    %
    % prospect theory weighting function
    %
    % a : scaling. increase = convex, < 0 = hyperbolic, 1 = linear
    % a > 0 & a < 1, averse to gambling (probs are lower than nominal)
    % a > 1, seeks gambling; probs are higher than nominal
    %
    % b : exponent.  increase: sigmoid shape, < 1, reverse slope
    % b > 0 & b < 1, overweight low prob, underweight high
    % b > 1, underweight low prob, overweight high

    fp = a .* (p.^b);
    wp = fp ./ (fp + (1-p).^b);

end

