function [X,names] = tor_get_param_mods(evtnums,rating,conditions,SPM)
% [X,names] = tor_get_param_mods(evtnums,rating,conditions,SPM)
%
% parametric modulators ACROSS sessions
% add call to this function after making design matrix
% then re-make design matrix...
% see model2_parammod in the rea dataset.
% works with SPM2
%
% enter event numbers you have param. modulators for, e.g., evtnums = [4 5 6]

% get param mod onsets

for k = evtnums
    parmods{k} = [];
    for j = 1:length(rating.run)
        parmods{k} = [parmods{k}; rating.run{j}.(conditions.names{k})'];
    end
end

for k = evtnums
    wh = find(isnan(parmods{k})); parmods{k}(wh) = nanmean(parmods{k});,

    % not necessary if we let SPM build the param mod.
    parmods{k} = parmods{k} - nanmean(parmods{k});
end

% set onsets to all onsets across sessions
SPMtmp = SPM;

%if we were going to create a stand-alone function
% put stuff into SPM
%SPM.Sess(1).row = 1:192;
%SPM.Sess(1).col = 1;
%SPM.Sess(1).Fc(1).i = 1;
%SPM.Sess(1).Fc(1).name = 'My_Regressor';


% cat onsets across conditions
for k = evtnums

    ons = [];
    for j = 1:length(SPMtmp.Sess),
        % adjust for this session
        offset = cumsum(SPMtmp.nscan(1:j-1));
        if isempty(offset),offset = 0;,else,offset=offset(end); ,end ,
        offset = offset * SPM.xY.RT;

        ons = [ons SPMtmp.Sess(j).U(k).ons + offset];,
    end

    % add onsets
    SPMtmp.Sess(1).U(k).ons = ons;

    % add param modulator
    SPMtmp.Sess(1).U(k).P.P = parmods{k};

    SPMtmp.Sess(1).U(k).P.name = 'other';
    SPMtmp.Sess(1).U(k).P.h = 1;

end


SPMtmp.nscan(1) = sum(SPM.nscan);


U = spm_get_ons(SPMtmp,1);

% use spm_get_bf
bf      = SPM.xBF.bf;

V   = 1;        % 1st-order Volterra.

[X,Xn,Fc] = spm_Volterra(U(evtnums),bf,V);

k   = SPMtmp.nscan(1);
fMRI_T     = SPM.xBF.T;
fMRI_T0    = SPM.xBF.T0;

% Resample regressors at acquisition times (32 bin offset)
%-------------------------------------------------------
try
    X = X([0:(k - 1)]*fMRI_T + fMRI_T0 + 32,:);
end

X = X(:,2:2:end); % get param mod cols only

% add X (param mod cols only) as user-specified regressors

indx = 1;
for k = evtnums
    names{indx} = {[U(k).name{1} '_PM']};
    indx = indx+1;
end


return

