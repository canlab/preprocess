function group_anova(params,A)
%
% by Tor Wager
%
% params    is a matrix of HRF estimates for all subjects
%           rows are estimates, trial types are concatenated in rows
%           columns are subjects, 1 col per subject
%
% A         is ANOVA structure from tor_setup_anova.m
%           contains indexes of which rows are which trial type, etc.

