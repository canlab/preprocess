
%% setup
x = [25 50 75 100 150 200 400 800 50 75 100 150 150 200 200]; x = [x x];
y = [0  0  0  0   0   0   0   0   25 50 50  50  100 100 150]; y = [y -y];
p = rand(size(x));

x = x';
y = y';
p = p';

truev = .7;
truea = .5;
trueb = .5;
truelo = 1.5;
%u = (x.^truev) .* wfcn(p,truea,trueb) + ((y .^ truev) .* (1-wfcn(p,truea,trueb))) + 5 .* randn(length(u),1);
%u is true u in head space.  put through inverse value function to go back
%to monetary value scale.

% true underlying utility
[trueu] = prospect_utility([truea trueb truev truelo],x,y,p);
% true report
ureport = vinv(trueu,truev,truelo);

 % x = best outcome
    % y = alternative outcome
    % p = prob. of getting x
    
    % a = scale, weighting function
    % b = exp., weighting function
    % m = exp., value function
    % lo = loss aversion
    
%% fit    
tor_fig(1,2)
[best_params,fit_u] = prospect_fit_data(x,y,p,ureport,1);
    
best_params
[truea trueb truev truelo]
[trueu] = prospect_utility([truea trueb truev truelo],x,y,p); [fit_u] = prospect_utility(best_params,x,y,p);
figure; plot(fit_u,trueu,'ko')

%% boot
% bootstrap to get variance
disp('Bootstrapping')
BOOT = bootstrp(50,@prospect_fit_data,x,y,p,u);
   
for i = 1:4
    subplot(1,4,i), hist(BOOT(:,i),50);
end

ste_params = std(BOOT);

%%
