function [meanest,min95est,max95est,ideal_data,meanfits,INFO,allb] = glm_estimate_dev(DX,snr,hrf,varargin)
% function [meanest,min95est,max95est,ideal_data,meanfits,INFO,allb] = glm_estimate_dev(DX,snr,hrf,[opt] ttype, niterations)
%
% DX is the model matrix containing a shape-free (time-domain deconvolution) model
% snr is the signal to noise ratio, defined here as resp height / noise variance
% hrf contains cell array of hrfs for each condition or single hrf
% [optional] ttype is which trial type to collect stats on
% 
%
%
    % rmsd      sqrt of mean squared dev. of estimates from true response
    % msdstd    single obs. 95% confidence interval for error (2-tailed)
    % msdb      abs dev. from ideal response for each time point
    % biasmean  mean bias (above or below true response) for each time pt.
    % meanest   mean estimate of the response (estimated hrf)
    % min95est  95% of individual regressions fall within this window...
    % max95est  between min95est and max95est
    %
    % assumes an ideal signal height of 1, so noise var = 1 / snr
    %
    % by Tor Wager
    
% optional argument is which trial type to test   
ttype = 1;
niterations = 1000;

if length(varargin) > 0
    if length(varargin) > 1
        niterations = varargin{2};
    end
    ttype = varargin{1};
end

[dummy, wd, wb] = DX_find_delta(DX);
mywin = [wd(ttype) : wd(ttype)+wb(ttype)];

% noise autocorrelation
x =1:100; y = 1 ./ (x + 7); y = y ./ max(y);
xc = y;

% for white noise:
%xc = [];


% for filtering
S = [];
%S= use_spm_filter(1.5,size(DX,1),'none','specify',64);
if ~isempty(S), DX = S * DX;,end


% set up hrf, if single function rather than cell array entered
if ~iscell(hrf),
    warning('Using single canonical HRF to fit - not convolved.')
    myh = hrf; clear hrf; for i = 1:length(wd), hrf{i} = myh;, end
end

%ideal_b = hrf{ttype};
%ideal_b = repmat(max(ideal_b),size(DX,2),1);
    
INFO.conv = 'all responses convolved with ideal';
INFO.condition_tested = ['Condition ' num2str(ttype)];
INFO.estlength = length(mywin);
INFO.num_condition_types = size(DX,2);
INFO.onsets = wd;
INFO.spacing = wb;
if isempty(xc), INFO.noisetype = 'white';, else, INFO.noisetype = 'colored';,end
INFO.autocorrelation = xc;
INFO.niterations = niterations;

% -------------------------------------------------------------------
% * set up multiple response
% -------------------------------------------------------------------

%ideal_data = [];
%for i = 1:size(ideal_delta,2)
%    mydata = conv(hrf{i},ideal_delta(:,i));
%    mydata = mydata(1:size(DX,1));
%    ideal_data(:,i) = mydata;
%end
ideal_data = sum(DX,2);
%ideal_data = ideal_data ./ var(ideal_data);

fprintf(1,'Iterating with snr %3.2f, saving trial type %3.0f, %3.0f iterations ',snr,ttype,niterations)
t1 = clock;

for i = 1:niterations  % iterations to generate population of estimates
  
    % -------------------------------------------------------------------
    % * set up noisy response
    % -------------------------------------------------------------------
    n = noisevector(size(DX,1),xc,var(ideal_data) ./ snr)';
    if ~isempty(S), n = S * n;, end
    noisy_ideal_data = ideal_data + n;
    b = pinv(DX) * noisy_ideal_data;

    %msdb(i,:) = (b - ideal_b)';
    fits(i,:) = (DX * b)';
    allb(i,:) = b';
    
    if mod(i,100) == 0, fprintf(1,'.'),end
end

fprintf(1,'%3.0f s elapsed\n',etime(clock,t1))

%rmsd = (mean((msdb .^ 2),2)).^.5;
%rmsd = sort(rmsd);
%msdstd(1,1) = rmsd(round(.025 * length(rmsd)));
%msdstd(1,2) = rmsd(round(.975 * length(rmsd)));
%rmsd = mean(rmsd);
    
%biasmean = mean(msdb);
%msdb = mean(abs(msdb));
meanest = mean(allb);
meanfits = mean(fits);
sallb = sort(allb);
min95est(1,:) = sallb(round(.025 * size(sallb,1)),:);
max95est(1,:) = sallb(round(.975 * size(sallb,1)),:);

return
