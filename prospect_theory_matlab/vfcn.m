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