function [pc,stats] = pca_npm2(X,perms,varargin);
% [pc,stats] = pca_npm2(X,perms,[MCD robust outlier removal],[input cov mtx not data])
%
% X     data matrix, columns are variables, rows observations
% MCD   optional: 1/0, do robust MCD trimming on data.  Default 0
% cov   optional: 1/0, flag for cov mtx input rather than data mtx.
% Default 0.
%
% stats.thresh = alphaval% level of eigvals, default 99%
% stats.wthresh = alphaval% level of eigval*w, divided back by w of actual pca
%
% example: 2-PC solution simulation
% X = rand(50,1); 
%X(:,2) = X(:,1) + rand(size(X,1),1);
%X(:,3) = X(:,1) + rand(size(X,1),1);
%X(:,4) = rand(size(X,1),1);
%X(:,5) = X(:,4) + rand(size(X,1),1);
%X(:,6) = X(:,4) + rand(size(X,1),1);
%[pc,stats] = pca_npm2(X,100);
%
% This version uses a step-down test to assess sig of each PC, then removes
% it and re-performs permutation test. 
% Works well for test example, but not for large numbers of similar
% variables (too many sig dimensions)

alphaval = 99;      % percentile cutoff for significance of eigenvalues

docov = 0;
if length(varargin) > 1, docov = varargin{2};,end

% ----------------------------------------------------------
% outliers
% ----------------------------------------------------------
dorob = 0;
if length(varargin) > 0, dorob = varargin{1};,end

if dorob
        [res]=fastmcd_noplot(X);
        
        % remove n most extreme outliers and recompute correlation
        wh = res.flag==0; nout = sum(res.flag==0);
        X(wh,:) = [];
        stats.nout = nout;
end

% ----------------------------------------------------------
% principal components
% ----------------------------------------------------------

if docov     % DATA IS COV MATRIX
% -----------------------------------------------------------------------
% permutations on covariance matrix -- don't know if this is right??
% -----------------------------------------------------------------------

    if any(diag(X)-mean(diag(X)) > eps),  warning('Thresholding may not be correct if diagonals are not all equal!'), end

    % WE should really permute pcares residuals for each dimension to get
    % the expected chance eigenvalues for each dimension
    
    [pc,eval,expl] = pcacov(X);
    score = NaN .* zeros(size(pc));
    eval = eval';
    
    Xvals = X - tril(X);   % permute these
    dia = diag(X); dia = diag(dia);

    for i = 1:perms,  
        [tmp,ind] = sort(rand(1,length(Xvals(Xvals~=0))));  % random index of nonzero elements in Xvals
        Xvtmp = Xvals(Xvals~=0);
        Xvtmp = Xvtmp(ind);                                 % random permutation of Xvals
        Xvt = Xvals; Xvt(Xvt ~= 0) = Xvtmp;                 % put randomized off-diagonal utl back into place
        
        Xvt = Xvt + Xvt' + dia;                              % construct full matrix
                                                            % diagonals
                                                            % maybe should
                                                            % be resorted
                                                            % as well but
                                                            % are not now.
        [pct,tmp] = pcacov(Xvt);
        cctst(i,:) = tmp';
    
        % for weights
        wtest(i,:) = max(pct * diag(tmp));
    end
        
    
else    % INPUT IS DATA
% -----------------------------------------------------------------------
% permutations on original data
% -----------------------------------------------------------------------
    
    % real-data solution

    [pc,eval,expl] = pcacov(corrcoef(X)); eval = eval'; expl = expl';
    [pc,score]=princomp(scale(X));
    pctst(:,1) = pc(:,1); evaltst = eval(1); expltst = expl(1);
    
    % for each PC, get residuals for PCs n - 1 and recompute
    % eigenvectors/values with others removed
    for i = 1:size(X,2)-1
        xtst = X;
        [res,rec] = pcares(xtst,i);
        [pctst1,evaltst1,expltst1] = pcacov(corrcoef(res));
        pctst(:,i+1) = pctst1(:,1); evaltst(i+1) = evaltst1(1); expltst(i+1) = expltst1(1);
    
    end
    
    stats.alphaval = (100 - alphaval) ./ 100;
    stats.evaltst = evaltst;
        
    %pc = -pc; score = -score;

    issig = 1;  % flag for significance of test -- to keep going or stop
    dimctr = 0; % dimension counter
    while issig
        
        % do permutations
        for i = 1:perms,   % permute columns and test H0: no relation
    
            xtst = X;
            % reshuffle 1st, then pcares--else pcares introduces
            % correlations
            for j = 1:size(res,2)       % reshuffle 
                xtst(:,j) = getRandom(xtst(:,j));
            end
            
            if dimctr > 0, [res,rec] = pcares(xtst,dimctr);, else, res = xtst;, end    % PCA residuals -- or X, for 1st pass


            [pct,tmp] = pcacov(corrcoef(res));
            cctst(i,1) = tmp(1);                % 1st eigenvalue of remaining data -- variances of PCs
    
            % for weights
            wtest(i,:) = max(pct * diag(tmp));

        end
    
        % test significance
        cctst = cctst(:,1);    % save 1st eigenvalues of permuted--for Ho distribution
        stats.thresh(dimctr + 1) = prctile(cctst,alphaval);
        if dimctr + 1 > length(evaltst), warning('All eigenvals sig! Something is wrong.');, keyboard, end
        
        issig = evaltst(dimctr + 1) > stats.thresh(dimctr + 1);  % test for significance. keep going if yes
        
        stats.sig(dimctr + 1) = issig;
        stats.p(dimctr + 1) = sum(cctst > repmat(evaltst(dimctr + 1),size(cctst,1),1)) ./ size(cctst,1);
        
        dimctr = dimctr + 1;
        
        
    end     % is-significant
   
end % if docov or dodata


% -----------------------------------------------------------------------
% save output and stuff
% -----------------------------------------------------------------------

stats.pc = pc;
stats.eigval = eval;
stats.expl = expl;

if perms
    %stats.cctst = cctst;
    stats.perms = perms;
    
 
    %t = 2; %input('Enter threshold eigenvalue');
    %stats.thresh = repmat(t,1,length(eval));
 
    stats.wh = find((1:length(stats.sig)) - cumsum(stats.sig) == 0);
    
       
    figure;plot(eval,'bo-','LineWidth',2); %hold on; plot(stats.thresh,'k.-')
    hold on; plot(eval(stats.wh),'ro-','LineWidth',2);
    
    legend({'Eigenvalues' 'Expected Ho eigenvalues'})
    %set(gca,'XTick',1:max(eval)),grid on
    
    
    stats.score = score;
    score = score(:,stats.wh);
    
    
    for i = 1:length(stats.wh),
        for j = 1:size(X,2)
            tmp = corrcoef(X(:,j),score(:,i));
            stats.wcor(j,i) = tmp(1,2);
        end
    end

    %stats.wthresh = prctile(wtest,95) ./ eval;
    stats.wthresh = repmat(.4,1,length(eval));
    pc = pc(:,stats.wh);
    
    % average X where weights are high or low on each
    class = []; classdata = [];
    for i = 1:length(stats.wh)
        tmp = stats.wcor(:,i) > stats.wthresh(i);
        if any(tmp)
            class(:,end+1) = tmp;
            classdata(:,end+1) = mean(X(:,find(class(:,end))),2);
        end
        
        tmp = stats.wcor(:,i) < -stats.wthresh(i);
        if any(tmp)
            class(:,end+1) = tmp;
            classdata(:,end+1) = mean(X(:,find(class(:,end))),2);
        end
    end
    
    stats.class = class;
    stats.classdata = classdata;
        
end

return
