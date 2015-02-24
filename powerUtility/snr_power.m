function [p,ps] = snr_power(fig)
% [p,ps] = snr_power(fig)
%
% the idea:
%
% r is the number of tests of a single voxel over n subjects.
% signal is added to r noise vectors of length n, 
% a t-test is performed on each of r vectors,
% correcting for v hypothetical voxels in the brain.
% (although only one hypothetical 'active' voxel is tested r times)
% and the probability of finding a positive result
% is computed as the number of positive results / r.
%
% false positive rates are computed using a signal strength of 0,
% and testing v voxels.  this will take a lot longer, i think.
%
% the same r noise vectors are used for various SNRs
% by adding different amounts of signal to it.
%
% this analysis is compared to the split analysis, which
% divides participants into two groups.
% the first group is used to identify ROIs at a nominal threshold,
% and the 2nd group is used to test those selected voxels,
% using a bonferroni correction based on the restricted space.
%
% power of each method is plotted on the same graph.
%
% Tor Wager, 2/14/02

% ------------------------------------------------------------------
% * primary variables
% ------------------------------------------------------------------
r   = 10000;        % number of tests to determine population
v   = 30000;        % no. of voxels in brain - multiple comparisons
t   = 1000;         % estimated number of truly active voxels
n   = 20;           % no. of subjects
snr = .5:.5:3;         % signal to noise ratio values to test
a1  = [.01 .001 .0001 .00001];          % nominal alpha value for 1st group in split test
verbose = 1;        % verbose output
dofull = 1;         % do all-data analysis
dosplit = 1;        % do split-data analysis
lcol = {'b--' 'g--' 'r--' 'c--'};        % line color for split analysis

if verbose, t1 = clock;, fprintf(1,'snr_power | setup | '),end

% ------------------------------------------------------------------
% * secondary variables
% ------------------------------------------------------------------
df      = n-1;                          % df for all-data analysis
a       = 0.05 / v;                     % alpha value for all-data an.
%tc      = tinv(1-(.05/v),df);           % critical t for all-data an.
pop     = randn(n,r);                   % normal noise


ns      = [n-floor(n/2) n-ceil(n/2)];   % no. of subjects in split test
dfs     = ns - 1;                       % df for split analysis
tc1     = tinv(1-(.05),dfs(1));         % crit t for ROI ident. in split
pop1    = pop(1:ns(1),:);               % first half of random data
pop2    = pop(ns(1)+1:end,:);           % 2nd half

if size(pop2,1) ~= ns(2), error('Pop computation error.'), end


% ------------------------------------------------------------------
% * testing all-data analysis over SNR
% ------------------------------------------------------------------

if dofull

    if verbose, t2 = clock;, fprintf(1,'all-data '),end

    p = test_pop(pop,snr,a,verbose);                % returns a vector of power at each snr val.

    if verbose, t3 = clock;, fprintf(1,'%3.0f s | ', etime(t3,t2)),end

end


% ------------------------------------------------------------------
% * testing split analysis over SNR
% ------------------------------------------------------------------
        
if dosplit
    for my_a = 1:length(a1)     % repeat analysis for each a1 threshold
    

        if verbose, t2 = clock;, fprintf(1,'split '),end

        indx = 1;
        for SNR = snr
            p1 = test_pop(pop1,SNR,a1(my_a),verbose);   % test of 1st group
            av       = p1 * t;                          % expected activated voxels in 1st group
            a2       = .05 ./ av;                       % bonferroni corrected critical p for 2nd group
                                                % based on the size of the active area in the 1st.
                                      
            ps{my_a}(indx) = test_pop(pop2,SNR,a2,verbose);   % test 2nd group with new alpha value
																							  % returns p of finding a voxel that survived 1st analysis
	    ps{my_a}(indx) = ps{my_a}(indx) .* p1;					   % returns p of finding a voxel in both analyses
            indx = indx + 1;
        end

        if verbose, t3 = clock;, fprintf(1,'%3.0f s |', etime(t3,t2)),end

    end
end

% ------------------------------------------------------------------
% * plot these results
% ------------------------------------------------------------------
if isempty(fig),fig = figure;, end
figure(fig)
hold on; set(gcf,'Color','w'), grid on

if dofull, plot(snr,p,'k','LineWidth',2), end
if dosplit, 
    for my_a = 1:length(a1)
        plot(snr,ps{my_a},lcol{my_a},'LineWidth',2), 
    end
end

ylabel('Power (prob. of finding positive result)','FontSize',14)
xlabel('Signal to noise ratio (mean / std)','FontSize',14)
title(['Power, n = ' num2str(n) ', ' num2str(t) ' truly active voxels'],'FontSize',14)

if verbose, t4 = clock;, fprintf(1,'\n done in %3.0f s total.', etime(t4,t1)),end

return




% ------------------------------------------------------------------
% * sub-functions
% ------------------------------------------------------------------

function [p,spop,cnt] = test_pop(pop,snr,a,verbose)
%
% assumes pop is null distribution, normal noise, u = 0, sd = 1.

r = size(pop,2);
index = 1;

for SNR = snr
    if verbose, fprintf(1,'.'),end
    testpop = pop + SNR;
    for i = 1:r
        h(i) = ttest((testpop(:,i)),0,a);       % test each group against mean 0, 2-tail
    end
    p(index) = sum(h) / r;                      % proportion of significant voxels
    index = index + 1;
    
end

                                                % for getting significant vectors only
                                                % don't really need this.
if nargout > 1                                  % use this with only 1 SNR!
    cnt = sum(h);
    if nargout > 2
        spop = pop(:,h);
        
    end
end

return
