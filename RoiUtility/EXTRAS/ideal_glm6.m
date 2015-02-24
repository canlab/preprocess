function [meanest,min95est,max95est,ALLINFO,hrf,snr,TR] = ideal_glm6(conditions,mspec,ttype)
% [meanest,min95est,max95est,ALLINFO,hrf,snr,TR] = ideal_glm6(conditions,mspec,ttype)
%
% tests glm model matrix directly against idealized data
% you put in the exact temporal sequence to be deconvolved,
% in the form of the X matrix.
%
% Tor Wager, 4/19/02
%
% inputs:
%	X	    glm design matrix
%	tp	    time points estimated for each condition in DX
%	TR	    repetition time of scan
%	ttype   trial types to test (out of 1:n different conditions in DX)
%           recommended for time saving to use ttype = a single number only
%
%
% this function is like ideal_deconv5, but tests variability across designs as well.

hrf = spm_hrf(mspec.TR);
hrf = hrf ./ max(hrf);

% hrfin contains cell array of hrfs for each condition
[X] = construct_model(mspec,conditions,[]);

[delta, wd, wb] = DX_find_delta(X);
for j = 1:length(wd)
    hrfin{j} = hrf;
    
    % when epoch or block regressor - now only works for one-part designs
    if length(conditions) >= j
        if isfield(conditions(j),'stimlength')
            if conditions(j).stimlength > 1
                hrfin{j} = conv(hrf,ones(conditions(j).stimlength,1));
                hrfin{j} = hrfin{j} ./ max(hrfin{j});
            end
        end
    end
    
end

% for custom hrf for one trial type
%hrf2 = sum([[zeros(5,1);hrf] [hrf;zeros(5,1)]],2);
%hrfin{2} = hrf2;


ndesigns = 100;
nnoise = 100;
snr = 1;
designs = 1:ndesigns;

for j = ttype
    
    clear rmsd, clear msdstd, clear msdb, clear biasmean, clear meanest, clear min95est, clear max95est
    
    for i = 1:ndesigns   % for this many designs
    
    % rmsd      sqrt of mean squared dev. of estimates from true response
    % msdstd    single obs. 95% confidence interval for error (2-tailed)
    % msdb      abs dev. from ideal response for each time point
    % biasmean  mean bias (above or below true response) for each time pt.
    % meanest   mean estimate of the response (estimated hrf)
    % min95est  95% of individual regressions fall within this window...
    % max95est  between min95est and max95est

        [X] = construct_model(mspec,conditions,[]);
        %ideal_delta = zeros(size(X,1),1);  
        %ideal_delta(round(conditions(1).onsets)+1) = 1;
        
        [meanest(i,:),min95est(i,:),max95est(i,:),ideal_data,meanfits(i,:),INFO] ...
        = glm_estimate_dev(X,snr,hrfin,j,nnoise);
    
        INFO.ttype = j;
    end

    
    
    % plotting
    
    
    
    plot_ideal_glm(meanest,min95est,max95est,ideal_data,meanfits,INFO,hrfin{j},designs,mspec.TR);

    
end

return
