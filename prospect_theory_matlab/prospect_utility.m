function [u] = prospect_utility(params,x,y,p)
    % [u] = prospect_utility(params,x,y,p)
    %
    % x = best outcome
    % y = alternative outcome
    % p = prob. of getting x
    % a = scale, weighting function
    % b = exp., weighting function
    % m = exp., value function
    %
    % Optional:
    % lo = loss aversion parameter, slope of negative value function
    %
    %

    if length(x) > size(x,1), x = x'; end
    if length(y) > size(y,1), y = y'; end
    if length(p) > size(p,1), p = p'; end
    
    a = params(1);
    b = params(2);
    m = params(3);
    lo = params(4);


    % x, y = two nominal values (Amounts)
    % p = nominal probability of getting it
    %
    % a, b = prob weighting params
    % m = value weighting param
    %
    % m = .7 , a = .5, b = .5

    % prob(x) * value(x) + prob(not_x) * value(y)
    % x is outcome, y is alternative outcome for choice

    u = wfcn(p,a,b) .* vfcn(x,m,lo) + ( (1 - wfcn(p,a,b)) .* vfcn(y,m,lo) );


end





function wv = vfcn(x,m,lo)
    % vfcn = inline('x.^m','x','m');
    % v = vfcn(x,m)
    % concave when x and m are both positive
    % creates loss of u relative to x when < 0

    if nargin < 3 || isempty(lo)
        % no loss aversion parameter
        lo = 1;
    end

    for i = 1:length(x)
        if x(i) < 0
            wv(i,1) = -lo .* abs(x(i)).^m;
        else
            wv(i,1) = x(i) .^ m;
        end
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
