function [rmsd,msdstd,msdb,biasmean,meanest,min95est,max95est,INFO] = dx_estimate_dev(DX,snr,hrf,varargin)
% function [rmsd,msdstd,msdb,biasmean,INFO] = dx_estimate_dev(DX,snr,hrf,[opt] ttype, [niterations],[acf], [truer])
%
% DX is the model matrix containing a shape-free (time-domain deconvolution) model
% snr is the signal to noise ratio, defined here as resp height / noise variance
% hrf contains cell array of hrfs for each condition or single hrf
% [optional] ttype is which trial type to collect stats on
% [optional] niterations is number of random noise vectors to test and average
% [optional] acf is autocorrelation function (vector of cross-correlation values)
%            empty = white noise
% [optional] truer vector of true response magnitudes
%
% ttype can be a number, a list of index values, or a vector
% indicating which trial types to get HRF estimates for and plot
% RIGHT NOW: this works with only one ttype at a time; iterate with ideal_deconv5, e.g., for multiple
%
% truer is a vector with one element per trial type,
% indicating the strength of 'true' activation in each DX to be convolved with 
% idealized "true" responses.
% Conditions are defined as the onsets of new trial types modeled with an FIR in DX.
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

[delta, wd, wb] = DX_find_delta(DX);    
mywin = [wd(ttype) : wd(ttype)+wb(ttype)];    

% noise autocorrelation
if length(varargin) > 2
    xc = varargin{3};
else
    x =1:100; y = 1 ./ (x + 7); y = y ./ max(y);
    xc = y;
end

% for white noise:
%xc = [];

if length(varargin) > 3
    truer = varargin{4};
else
    truer = ones(1,size(delta,2));
end

% for filtering
S = [];
%S= use_spm_filter(1.5,size(DX,1),'none','specify',64);
if ~isempty(S), DX = S * DX;,end


% set up hrf, if single function rather than cell array entered
if ~iscell(hrf),
    myh = hrf; for i = 1:length(wd), hrf{i} = myh;, end
end

ideal_b = hrf{ttype}(1:length(mywin)) .* truer(ttype) ;

INFO.conv = 'Trial onsets convolved with canonical HRF; magnitudes in truer';
INFO.truer = truer;
INFO.condition_tested = ['Condition ' num2str(ttype)];
INFO.ttype = ttype;
INFO.hrf = hrf;
INFO.estlength = length(mywin);
INFO.num_condition_types = delta;
INFO.onsets = wd;
INFO.spacing = wb;
if isempty(xc), INFO.noisetype = 'white';, else, INFO.noisetype = 'colored';,end
INFO.autocorrelation = xc;
INFO.niterations = niterations;
INFO.delta = delta;

% -------------------------------------------------------------------
% * set up multiple response
% -------------------------------------------------------------------

ideal_data = [];
for i = 1:size(delta,2)
    mydata = conv(hrf{i},delta(:,i)) .* truer(i);
    mydata = mydata(1:size(DX,1));
    ideal_data(:,i) = mydata;
end
ideal_data = sum(ideal_data,2);

INFO.ideal_data = ideal_data;

fprintf(1,'Iterating with snr %3.2f, saving trial type %3.0f, %3.0f iterations ',snr,ttype,niterations)
t1 = clock;

for i = 1:niterations  % iterations to generate population of estimates
  
    % -------------------------------------------------------------------
    % * set up noisy response
    % -------------------------------------------------------------------
    n = noisevector(size(DX,1),xc,1 ./ snr)';
    if ~isempty(S), n = S * n;, end
    noisy_ideal_data = ideal_data + n;
    b = pinv(DX) * noisy_ideal_data;

    msdb(i,:) = (b(mywin) - ideal_b)';
    allb(i,:) = b(mywin)';
    
    if mod(i,100) == 0, fprintf(1,'.'),end
end

fprintf(1,'%3.0f s elapsed\n',etime(clock,t1))

rmsd = (mean((msdb .^ 2),2)).^.5;
rmsd = sort(rmsd);
msdstd(1,1) = rmsd(round(.025 * length(rmsd)));
msdstd(1,2) = rmsd(round(.975 * length(rmsd)));
rmsd = mean(rmsd);
    
biasmean = mean(msdb);
msdb = mean(abs(msdb));
meanest = mean(allb);

sallb = sort(allb);
min95est(1,:) = sallb(round(.025 * size(sallb,1)),:);
max95est(1,:) = sallb(round(.975 * size(sallb,1)),:);

return
