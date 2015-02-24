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

 
