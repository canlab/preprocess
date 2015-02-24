function v = uinv(u,m)
% inverse value function
%
% v(u) = subject report cash equivalence
% because subject reported utility is weighted by value function
% compare vu with subject reported utility
%vu = vinv(u,m,lo);
%
% RP = reported price =
% vinv(u) = vinv(prospect_utility(params, options))


whpos = u > 0;
whneg = ~whpos;
v = zeros(size(u));

v(whpos) = u(whpos) .^(1./m);

v(whneg) = -((-u(whneg)).^(1./m));




end