function [pc,stats] = pca_npm(X,perms,varargin)
% [pc,stats] = cancor(X,perms,[MCD robust outlier removal],[input cov mtx not data],[make figure])
%
% X     data matrix, columns are variables, rows observations
% MCD   optional: 1/0, do robust MCD trimming on data.  Default 0
% cov   optional: 1/0, flag for cov mtx input rather than data mtx.
% Default 0.
%
% stats.thresh = 95% level of eigvals
% stats.wthresh = 95% level of eigval*w, divided back by w of actual pca
%
% [pc,stats] = pca_npm(X,100,0,0,0); no robust, enter data matrix, no new
% figure
% tor wager

docov = 0; dofig = 1;
if length(varargin) > 1, docov = varargin{2};,end
if length(varargin) > 2, dofig = varargin{3};,end

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
        [pct,tmp] = pcacov_t(Xvt);
        cctst(i,:) = tmp';
    
        % for weights
        wtest(i,:) = max(pct * diag(tmp));
    end
        
    
else    % INPUT IS DATA
% -----------------------------------------------------------------------
% permutations on original data
% -----------------------------------------------------------------------
    
    [pc,eval,expl] = pcacov_t(corrcoef(X)); eval = eval'; expl = expl';
    [pc,score]=princomp_t(scale(X));
    %pc = -pc; score = -score;


    for i = 1:perms,   % permute columns and test H0: no relation
    %disp(['perm ',num2str(i)]);
        xtst=X;
        for j = 1:size(X,2)
            xtst(:,j) = getRandom(X(:,j));
        end
        [pct,tmp] = pcacov_t(corrcoef(xtst));
        cctst(i,:) = tmp';
    
        % for weights
        wtest(i,:) = max(pct * diag(tmp));

    end
   
end % if docov


% -----------------------------------------------------------------------
% save output and stuff
% -----------------------------------------------------------------------

stats.pc = pc;
stats.eigval = eval;
stats.expl = expl;

if perms
    %stats.cctst = cctst;
    stats.perms = perms;
    try
        stats.thresh = prctile(cctst,95);
    catch
        stats.thresh = prctile_t(cctst,95);
    end
    
    stats.sig = eval > stats.thresh;
    stats.p = sum(cctst > repmat(eval,size(cctst,1),1)) ./ size(cctst,1);
    stats.wh = find((1:length(stats.sig)) - cumsum(stats.sig) == 0);
    
    
    if dofig,figure;,end
    plot(eval,'ko:','LineWidth',1); 
    if dofig, hold on; plot(prctile(cctst,95),'k.-'), end
    hold on; plot(eval(stats.wh),'ko-','LineWidth',2,'MarkerFaceColor','r');
    legend({'Eigenvalues' 'Significant eigenvalues'})
    %set(gca,'XTick',1:max(eval)),grid on
    
    %t = 2; %input('Enter threshold eigenvalue');
    %stats.thresh = repmat(t,1,length(eval));
    

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
