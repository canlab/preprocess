%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Obligatory first part of the automation script                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  fmriTEST = 1 to run the simulation mode 
%           = 0 to run the actual post-processing 
%  fmriDIR  = the directory in which all your data is situated 
%----------------------------------------------------------------------- 
global fmriTEST fmriDIR;             %set global variabels 
fmriTEST = 1; 
fmriDIR  = '/images_spike1/beatse/sidev6'; 
f99_CWD('LOG');                   %make the LOG and RESULTS directory 
  

% calculation of all timing parameters 
% ====================================================================== 
onset{1} = [0 20 40 60 80 100 -1 -1 -1 -1  -1  -1;... 
            2  6 14 44 50  58 66 70 76 83  99 106;... 
            4 10 20 30 38  47 62 68 78 88 102 112]; 
onset{2} = [10 30 40 70 90 105 -1 -1 -1 -1  -1  -1;... 
             4  8 19 42 55  56 66 70 78 85 100 108;... 
             6 12 22 31 38  40 60 66 76 84 104 112]; 

varab{1} = [0 0 0 0 0 0 0 0 0 0 0 0;... 
             1  5  1  7  2  8  2  4  2  5  2  6;... 
            0 0 0 0 0 0 0 0 0 0 0 0]; 
varab{2} = [0 0 0 0 0 0 0 0 0 0 0 0;... 
             1  10  1  2  2  10  2  4  2  5  2  6;... 
            0 0 0 0 0 0 0 0 0 0 0 0]; 

length{1} = [10 0 0]; 
length{2} = [15 0 0]; 


% additional parameters 
% ====================================================================== 
names = ['sid1';'sid2';'sid3']; 
% user defined hrf function 
AA = [1:147]; 
BB = sin(7*pi\AA); 
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% function f99_stat' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% ====================================================================== 
% 'general' 
% ====================================================================== 
%  - repetition time 
%  - number of scans per session 
%  - are sessions replicated exactly (1=yes / 0=no) 
%  - number of conditions or trials 
%  - matrix with names of conditions/trials 
%  - directory where to write results 
% ---------------------------------------------------------------------- 
f99_stat_general(... 
    3.5,[105 105],0,... 
    3,names,... 
    'test') 
  

% ====================================================================== 
% 'stochastic design' 
% ====================================================================== 
%  - DO YOU WANT A STOCHASTIC DESIGN ('yes'/'no') 
%      * if you answer yes here, answer no in the next part!! 
%      * if you answer no here, no further input arguments are required! 
%  - include null event? (1=yes / 0=no ; for each session) 
%  - soa for each session 
%  - matrix with occurence probabilities for each trial type 
%     (nsess x # trial) 
%  - stationary (1) or modulated (0) occurence probabilities 
%     (for each session) 
% ---------------------------------------------------------------------- 
f99_stat_stoch(... 
    'no') 
  

% ====================================================================== 
% 'non stochastic design' 
% ====================================================================== 
%  - DO YOU WANT A DETERMINISTIC DESIGN ('yes'/'no') 
%      * answer yes here only if you said no in the previous part 
%      * if you answer no here, no further input arguments are required! 
%  - type of SOA ('Fixed' or 'variable') 
%  - time to first trial (scans) 
%         if 'fixed' : onset of first trial presentation 
%         if 'variable' : onset times of all trials 
%            * only one session or sessions are replicated 
%              onset = [n_trials*n_trial_types] 
%            * several sessions, not replicated 
%              onset{1}      = [n_trials*n_trial_types] 
%              onset{2}      = [n_trials*n_trial_types] 
%              ... 
%              onset{n_sess} = [n_trials*n_trial_types] 
%  - variable durations for all sessions and conditions 
%         (ncond x # onsets for each cond) for each session 
%  - soa for each trial type 
%         only if 'fixed' : 
% ---------------------------------------------------------------------- 
f99_stat_nonstoch(... 
    'yes',... 
    'Variable',... 
    onset,... 
    varab,... 
    []) 
  

% ====================================================================== 
% 'respons modulation' 
% ====================================================================== 
%  - parametric modulation? ('none','time','other') 
%  - name of parameter 
%         (in case of 'other') 
%  - expansion type ('linear','exponen','polynom') 
%  - time constant 
%         (in case of 'exponen') 
%  - order of polynomial expansion 
%         (in case of 'polynom') 
%  - which trials to modulate 
%  - matrix with parameters for each selected covariate 
%    [#starting points per covariate x #selected covariates] 
% ---------------------------------------------------------------------- 
f99_stat_modul(... 
    'none',... 
     '',... 
     'exponen',... 
     400,[],... 
     [1 2 3],... 
     []) 
  

% ====================================================================== 
% 'basisfunctions to use' 
% ====================================================================== 
%  - matrix with types of trials 
%      * event ('ev') or epoch ('ep') 
%      * one string per trial type (e.g. ['ep';'ep';'ev']) 
%  - user specified hrf function of this subject 
%      length of this function can be checked using fmriTEST 
%  - model type for covariates of interest 
%      1 = basis functions  (Discrete Cosine Set 
%      2 = basis functions  (Mean & exponential decay) 
%      3 = fixed response   (Half-sine) 
%      4 = fixed response   (Box-car) 
%  - model type for event-related responses 
%      5 = hrf (alone) 
%      6 = hrf (with time derivative) 
%      7 = hrf (with time and dispersion derivatives) 
%      8 = basis functions (Fourier set) 
%      9 = basis functions (Windowed Fourier set) 
%     10 = basis functions (Gamma functions) 
%     11 = basis functions (Gamma functions with 
%          derivatives) 
%     * e.g. [4 4 6] 
%  'event' 
%  - window length          (for Cov=8 Cov=9) 
%  - order of Fourier set   (for Cov=8 Cov=9) 
%  'epoch' 
%  - number of basis functions (if Cov=1) 
%  - convolve with HRF (1=yes / 0=no) 
%  - add temporal differences (1=yes / 0=no) 
%  - epoch length (scans) for each condition 
%  'additional' 
%  - interactions among trials (Volterra series)? (1=yes / 0=no) 
%  - matrix with user specified regressors (#regr x nscan) 
%  - matrix with strings of names of regr_us 
% ---------------------------------------------------------------------- 
f99_stat_basisfunct(... 
    ['ep';'ev';'ev'],... 
    [],... 
    [4 5 5],... 
    0,0,... 
    0,1,0,length,... 
    0,... 
    [],... 
    []) 

% ====================================================================== 
% 'construction of the design matrix' 
% ====================================================================== 
f99_make_DesMtx 
  

% ====================================================================== 
% 'estimation of the specified model' 
% ====================================================================== 
%  - [matrix of directories] 
%  - [matrix of experiment names] 
%  - directory where to find SPM_fMRIDesMtx 
%    directory where to write the results 
%  - global normalisation ('Scaling' or 'None') 
%  - burst mode? (1=yes, 0=no) 
%  - apply Hipgh Pass Filter? ('none' or 'specify' or 'auto') 
%  - if 'specify': cut off period for each session [1*nsess] 
%  - Low Pass Filter? ('none', 'Gaussian' or 'hrf') 
%  - user defined Gausssian FWHM (secs) 
%  - intrinsic correlations ('none' or 'AR(1)') 
%  - calculate trial specific F-contrasts (1=yes / 0=no) 
% ---------------------------------------------------------------------- 
f99_stat_calc(... 
    ['run01';'run02'],... 
    ['sn*.img';'sn*.img'],... 
    'test',... 
    'Scaling',0,... 
    'specify',[120 120],... 
    'hrf',[],... 
    'none',... 
    0) 


