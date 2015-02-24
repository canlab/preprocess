function [X,hrf] = get_X(EXPT,v,x,y,z,olen)
% [X,hrf] = get_X(EXPT,v,x(i),y(i),z(i),olen)
%
% Tor Wager
% Returns a model matrix X defined by estimated HRFs
% stored in matrix v (from dx_beta_imgs)
%
% Each predictor's HRF is the centered HRF for EACH
% condition.  Not appropriate for testing diffs between
% conditions.
% 
% if the HRF is the AVERAGE
% HRF across conditions, this is appropriate for testing
% for differences between conditions.
%
% 
%
% Called from get_hrf_X.m

% ----------------------------------------------------
% * get HRF and build model
% ----------------------------------------------------

for i = 1:length(EXPT.DX.dxtrialonsets), 
    hrf{i} = squeeze(v(x,y,z,EXPT.DX.dxtrialonsets(i):EXPT.DX.dxtrialonsets(i)+EXPT.DX.numframes-1));
    
    m = conv(hrf{i},EXPT.DX.DX{1}(:,EXPT.DX.dxtrialonsets(i)));
    X(:,i) = m(1:olen);
    X(:,i) = X(:,i) - mean(X(:,i));
    
end

X(:,end+1) = 1;

return