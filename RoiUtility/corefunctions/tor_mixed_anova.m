function [sig,T] = tor_mixed_anova(data,GROUP,factornames)
% tor_mixed_anova(data,GROUP,factornames)
%
% Tor Wager
% 
% Runs the following:
% [P,T,STATS,TERMS]=anovan(Y,GROUP,'full',3,gname,'off')
% [P,T,STATS,TERMS]=anovan(Y,GROUP,MODEL,SSTYPE,GNAME,DISPLAYOPT)
% Subject is fixed effect, unless you use this modified function
%
% my current understanding:
% since subject is a random effect, and the other factors are
% repeated measure fixed effects, the error term should be based on
% the subject * condition interaction (highest-order interaction)
% Running the anova without the subject term, the mse is effectively
% the subject * condition interaction, so we just have to adjust the
% dfe to reflect the repeated measures design.
%
% MSB / MSE = F, so replace MSE with sub * condition interaction
%
% data = matrix of individual subject parameter values
%   rows are observations (time points in each trial type, trial
%   types concatenated)
%   columns are subjects
%
% GROUP is as defined in help anovan
%   can be created with tor_setup_anova_ui.m
%   GROUP{1} = subject
%   (subject not included in anovan command)
%
% inputs, except for data, can be defined using tor_setup_anova_ui.m
% Example:
% [sig,T] = tor_mixed_anova(params{i},EXPT.ANOVA.factorcodes,EXPT.ANOVA.factornames);

y = reshape(data,size(data(:),1),1);

[P,T]=anovan(y,GROUP(2:end),'full',3,factornames(2:end),'off');

% adjust the error degrees of freedom because factors are w/i subject effects
dfe = T{end-2,3} * (max(GROUP{1}) - 1);


%dfe = T(end-2,3);   % dfe is the df for highest-order interaction
%mse = T(end-1,5);   % mse should be the same as before...sub * cond interaction

% recalculate F and p values using highest-order interaction
for i = 2:size(T,1) - 2
    T{i,7} = 1-fcdf(T{i,6},T{i,3},dfe);
    
    if T{i,7} < .05, sig(i-1) = 1;, else sig(i-1) = 0;,end
end

return