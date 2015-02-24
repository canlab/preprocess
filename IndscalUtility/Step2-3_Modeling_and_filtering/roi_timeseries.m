function o = roi_timeseries(m,o)
% o = roi_timeseries(m,o)
%
% RoI summary measures for voxels in m
% 
% For each column of m:
%
% - remove session means
% - high-pass filter
% - Windsorize to 3 std, recursive, with spike correction
% - Smooth (Filter) each session
% - Scale session variances to 1
%       - this makes sessions contribute equally to multivar measures
%         and to correlations with other timeseries.
%       - otherwise 1 session with high variance can affect covs
%       disproportionately
%
% Then compute summary measures:
% - average
% - significant principal components
% - factor scores
% - independent components (on all voxels, m, not reduced PCA space)
%
% Then deconvolve with canonical HRF
% or your estimate, based on input function
%
% You put in filtering options:
% 		o.scanadjust = 1;
%		o.percent = 0;
%		o.TR = .5;
%		o.nscans = length(RoI.ts.avg) ./ 860;
%        o.nruns = length(RoI.ts.avg) ./ 860;
%       o.filtertype = 'spm';
%	    o.HP = 82;			% only necessary for SPM
%	    o.doLP = 0;
%       o.doHP = 1;
%       o.trimts = 3;
%       o.docustomadjust = ***
% RECoMMENDED: no percent, no doLP      
%
% Varimax rotation tries to ensure that a subset of voxels load on a
% component, and each component is a combination of voxels that load highly
% or very little.
% Thus, voxels can load on more than one component, but each component is
% selected to come from only a subset of voxels.
% Component scores can be correlated to some degree.
%
% FastICA on rows (m', voxels x time) returns orthogonal time courses
% in fastICA, rows are 'signals'

dplot = 0;
fprintf(1,'roi_timeseries.setup')

% defaults

o.trimts = 3; % trim to n std
o.doLP = 0;
o.doHP = 1;
o.filtertype = 'spm';
o.percent = 0;
o.scanadjust = 1;

nperscan = size(m,1) ./ o.nscans;
o.y = m(:,1);

if isfield(o,'xpx')
    xpx = o.xpx;
    S = o.S;
    
else
    % do first one to set up filter, etc.
    [tmp, o, X, S] = filterAdjust(o);
    xpx = X * pinv(X);
    o.xpx = xpx;
    o.S = S;
    
end


%figure('Color','w'); subplot(2,1,1);plot(o.y); title('Original'),subplot(2,1,2); plot(tmp,'r'),title('Filtered'),drawnow



% process each voxel independently
fprintf(1,'.filter')

for i = 1:size(m,2)
    
    tmp = m(:,i);
    
    %tmp = (tmp - xpx * tmp); % remove session means - too slow!!
    
    % HP filter each session
    tmp = reshape(tmp,nperscan,o.nscans);
    tmp = tmp - repmat(mean(tmp),size(tmp,1),1);
    tmp = S * tmp;
    tmp = scale(tmp);       % make all sess have same variance 
    tmp = tmp(:);
    tmp = trimts(tmp,o.trimts,[],1);
    
    m(:,i) = tmp;
    
end

if dplot, figure;plot(m);, title('All voxels'),drawnow, end
% average, etc.

fprintf(1,'.Cond# = %3.1f.',cond(cov(m)));

o.average = nanmean(m')';
o.m = m;
fprintf(1,'.PCA')
[pcs,stats] = pca_npm(m,5);     % nonpar thresholding for significant pcs
o.nfact = sum(stats.sig,2);
o.pcscore = stats.score(:,1:o.nfact);
fprintf(1,'.FA')
try
    %lam = factoran(m, o.nfact);
    [lam, uniqueness, rotmtx, STATS, F] = factoran(m,o.nfact);
    o.varimax = F;
    %o.varimax = m * lam;   % gives diff estimates, altho they match ICA
    %better than the F ones...
    o.lam = lam;
catch
    fprintf(1,'No factor structure obtained.')
    o.varimax = NaN .* zeros(size(o.pcscore,1),o.nfact);
    o.lam = NaN;
end
   
fprintf(1,'.ICA')
try
    o.ica = fastICA(o.m','lastEig',o.nfact,'verbose','off','displayMode','on','maxNumIterations',100)';
catch
    fprintf(1,'(err)')
end

fprintf(1,'\n')

if dplot, 
    try
        cla, close, hold on; 
        figure('Color','w');
    for i = 1:o.nfact
         subplot(o.nfact,1,i); hold on;
        if i == 1
            plot(scale(o.average),'k'); 
        end
        plot(scale(o.pcscore(:,i)),'b'); plot(scale(o.varimax(:,i)),'g'); plot(scale(o.ica(:,i)),'r');
        if i == 1,legend({'Average' 'PCA' 'Varimax' 'ICA'}),else, legend({'PCA' 'Varimax' 'ICA'}),end
        title(['Scaled canonical timeseries: F' num2str(i)]);
    end
    
    catch
    end

end
    
return



%t = cov(m);
%[pc,latent,exp] = pcacov(t);
%[pc,score,lat] = princomp(m);
%o.pcscore = score(:,1:2);


% deconvolution

cf{1} = condf;
[DX,sf] = tor_make_deconv_mtx(cf,30,1);
DX(:,end+1) = 1;
%o.hrf = pinv(DX) * o.average;

S.data = o.average; S.condf = condf; S.color = {'r'}; S.method = 'average'; S.window = [0;30]; S.ste = 1; S =trialavg(S);
o.hrf = scale(S.avg{1});

C = (convmtx(o.hrf',nperscan));
C = C(1:size(C,2),:);
C = pinv(C);

tmp = o.average;
for j = 1:o.nscans
    tmp((j-1) .* nperscan + 1:j.*nperscan) = C * tmp((j-1) .* nperscan + 1:j.*nperscan);
end
o.dxaverage = tmp;

