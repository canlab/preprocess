function v = vinv(u,m,lo)
% inverse value function
%
% v(u) = subject report cash equivalence
% because subject reported utility is weighted by value function
% compare vu with subject reported utility
%vu = vinv(u,m,lo);
%
% RP = reported price =
% vinv(u) = vinv(prospect_utility(params, options))

if nargin < 3 || isempty(lo)
    % no loss aversion parameter
    lo = 1;
end

whpos = u > 0;
whneg = ~whpos;
v = zeros(size(u));

v(whpos) = u(whpos) .^(1./m);

v(whneg) = -((-u(whneg)./lo).^(1./m) );




end