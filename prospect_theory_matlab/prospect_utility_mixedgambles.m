function u = prospect_utility_mixedgambles(params,x,y,p)
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
    
    %for use with GA
%     a = params{1};
%     b = params{2};
%     m = params{3};
%     lo = params{4};
    
    u = zeros(length(x),1);
    wvalue = zeros(length(x),1);

    % x, y = two nominal values (Amounts)
    % p = nominal probability of getting it
    %
    % a, b = prob weighting params
    % m = value weighting param
    %
    % m = .7 , a = .5, b = .5

    % prob(x) * value(x) + prob(not_x) * value(y)
    % x is outcome, y is alternative outcome for choice
    for i = 1:length(x),
        w_xval(i,1) = vfcn(x(i),m,lo);
        w_yval(i,1) = vfcn(y(i),m,lo);
        wpm = wfcn(1-p(i),a,b); %for mixed gamble
        wp = wfcn(p(i),a,b); %for nonmixed gamble
        
        % equation 3
        if x(i).*y(i) <0,   % if mixed gamble
          u(i) = w_xval(i).*wp + w_yval(i).*wpm;
        else % not a mixed gamble 
          u(i) = w_xval(i).*wp + w_yval(i).*(1-wp);
        end
    end

  


end




% equation 1
function wv = vfcn(x,m,lo)
    % vfcn = inline('x.^m','x','m');
    % v = vfcn(x,m)
    % concave when x and m are both positive
    % creates loss of u relative to x when < 0

        if x < 0
            wv = -lo.* abs(x).^m;
        else
            wv = x .^ m;
        end
       

end



% equation 2
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
