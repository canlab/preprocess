%mvroi_generate_design_plugin

%% generate design matrices %%%%%%%%%%%%%%%%
%
% Needs:
% DATA.onsets           onsets of event times in SPM-style format (cell array
%                       of times, condition nested within session
% numsub                number of subjects
% DATA.SPEC.firpoints   vector or integer of number of betas to estimate
%                       for each condition; same as EXPT.DX.numframes
% DATA.SPEC.spersess    Number of images in each session, [n1 n2 ... nx]
%
%
% Will take:
% DATA.DX        supersedes DATA.onsets
%
% Creates:
% DATA.DX        cell array of one design matrix per subject
%
% Uses functions (nested to reflect structure of calls):
% build_fir_model
%   onsets2delta
%   tor_make_deconv_mtx3
%       intercept_model
%       
numsub=DATA.SPEC.numsub;

% look for existing design matrices
if isfield(DATA,'DX'),
    fprintf(1,'Found existing design matrices in DX . ');
    dobuild = input('Re-create from onsets? (1/0) : ');
else
    dobuild = 1;
end

if ~dobuild
    % already got em, don't rebuild
    if length(DATA.DX) < numsub
        fprintf(1,'*** DX seems to short-replicating first matrix for each subject . \n');
        for s = 1:numsub, DATA.DX{s} = DATA.DX{1};, end
    end
 
else
 % build design matrices  
 
    if ~isfield(DATA.SPEC,'onsets')
        error('Enter onset vectors in DATA.SPEC.onsets or use already-created DATA.DX cells');
    end
     
    % check for length of trigs (onsets)
    if length(DATA.SPEC.onsets) < numsub
        fprintf(1,'*** Fewer onset vector cells than subjects - replicating onsets for each subject . \n');
        tmp = DATA.SPEC.onsets;
        for s = 1:numsub, DATA.SPEC.onsets{s} = tmp;, end
        clear tmp
    end
        
    fprintf(1,'generating design matrices . ');
    for s=1:numsub                                  % for all subjects
        DATA.DX{s} = build_fir_model(DATA.SPEC.onsets{s},DATA.SPEC.firpoints,DATA.SPEC.spersess);     
    end
    
end
