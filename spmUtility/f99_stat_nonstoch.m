function f99_stat_nonstoch(N_STOC,Sstr,onset,durat,soa_n_stoc)

% 'non stochastic design'
% DO YOU WANT A DETERMINISTIC DESIGN ('yes'/'no')
%          answer yes here only if you said no in the previous part
% Sstr        - type of SOA ('fixed' or 'variable')
% onset       - time to first trial (scans)
%               if 'Fixed' : onset of first trial presentation
%               if 'Variable' : onset times of all trials
%                  * only one session or sessions are replicated
%                    onset = [n_trials*n_trial_types]
%                  * several sessions, not replicated
%                    onset{1}      = [n_trials*n_trial_types]
%                    onset{2}      = [n_trials*n_trial_types]
%                    ...
%                    onset{n_sess} = [n_trials*n_trial_types]
% durat       - variable durations for all sessions and conditions 
%                   (ncond x # onsets for each cond) for each session 
% soa_n_stoc  - soa for each trial type
%               only if 'Fixed' : 

global fmriTEST fmriDIR kulRES
logfile=[fmriDIR filesep 'LOG' filesep 'design.log'];

if strcmp(deblank(lower(N_STOC)),'no')
    Sstr       = '';
    onset      = [];
    soa_n_stoc = [];
else
    f99_log(logfile,['DETERMINISTIC']);
    f99_log(logfile,['  SOAs are ' Sstr]);
    if ~isempty(durat)
      f99_log(logfile,['  epochs with variable durations']);
    end
end

filename = [fmriDIR filesep 'RESULTS' filesep kulRES filesep 'KUL_stat_nonstoch'];
str = ['save ' filename ' Sstr onset durat soa_n_stoc'];
eval(str)
