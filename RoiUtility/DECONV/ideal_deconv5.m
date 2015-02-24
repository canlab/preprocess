function [ALLINFO,rmsd,msdstd,msdb,biasmean,meanest,min95est,max95est,hrf,snr,TR] = ideal_deconv5(DX,TR,ttype,varargin)
% [INFO,rmsd,msdstd,msdb,biasmean,meanest,min95est,max95est,hrf,snr,TR] = ideal_deconv5(DX,TR,ttype,[acf],[truer])
%
% tests deconvolution matrix directly against idealized data
% you put in the exact temporal sequence to be deconvolved,
% in the form of the DX matrix.
%
% Tor Wager, 4/19/02
%
% inputs:
%	DX	    deconvolution matrix
%	tp	    time points estimated for each condition in DX
%	TR	    repetition time of scan
%	ttype   trial types to test/plot (out of 1:n different conditions in DX)
%           ttype can be a number, a list of index values, or a binary logic vector
%           
%           Conditions are defined as the onsets of new trial types modeled with an FIR in DX.
%   acf     Autocorrelation function vector [optional]
% % truer   vector of true response magnitudes [optional] 
%
% Output: Bias, error, deconv plots

if length(varargin) > 0, acf = varargin{1};, else, acf = [];, end

hrf = spm_hrf(TR);
hrf = hrf ./ max(hrf);

% hrfin contains cell array of hrfs for each condition
[delta, wd, wb] = DX_find_delta(DX);
for j = 1:length(wd)
    hrfin{j} = hrf;
end

if length(varargin) > 1, truer = varargin{2};, else, truer = ones(1,size(delta,2));, end

% for custom hrf for one trial type
%hrf2 = sum([[zeros(5,1);hrf] [hrf;zeros(5,1)]],2);
%hrfin{2} = hrf2;

snr = [.3 .4 .6 .9 1.3 1.8];
snr = [.5 .9 1.1 1.5 1.9];

for j = ttype
    
    clear rmsd, clear msdstd, clear msdb, clear biasmean, clear meanest, clear min95est, clear max95est
    
    for i = 1:length(snr)  
    
        % rmsd      sqrt mean squared dev. of estimates from true response
        % msdstd    single obs. 95% confidence interval for error
        % msdb      abs dev. from ideal response for each time point
        % biasmean  mean bias (above or below true response) for each time pt.

        [rmsd(i),msdstd(i,:),msdb(i,:),biasmean(i,:),meanest(i,:),min95est(i,:),max95est(i,:),INFO] ...
        = dx_estimate_dev(DX,snr(i),hrfin,j,200,acf,truer);
    
    end

    plot_ideal_deconv5(rmsd,msdstd,msdb,biasmean,meanest,min95est,max95est,INFO,hrfin{j},snr,TR,truer);

    INFO.TR = TR;
    ALLINFO(j) = INFO;

end

return
