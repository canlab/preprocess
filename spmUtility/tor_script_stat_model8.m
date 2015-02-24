%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Obligatory first part of the automation script                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  fmriTEST = 1 to run the simulation mode 
%           = 0 to run the actual post-processing 
%  fmriDIR  = the directory in which all your data is situated 
%----------------------------------------------------------------------- 
global fmriTEST fmriDIR;             %set global variabels 
fmriTEST = 0; 
fmriDIR  = '/data1/intext';
f99_CWD('LOG');                   %make the LOG and RESULTS directory 

% model2: 6 regressors
% model3: volterra
% model5: left/right
% model6: 6 regressors, first 4 scans only
% model8 : SPM regressors, 6 per scan
% model9: SPM regs, myVi, no temporal smoothing
% model 10: SPM regs, 6 per, first half
% model 11: SPM regs, 6 per, 2nd half
% model 12: same as 8, but adding movement parameters.
% model 13: same as 12, but with temporal derivative added

disp('============================ tor setup ==================================')
dostat = 0;
docontrasts = 1;
donorm = 1;
dorfxcopy = 0;
checkindss = 0;
runfixedfx = 0;
model = '8';
numscans = 8;
subs = [9];
ss = subs; modelnum = model;
drive = [2];


if dostat
   
disp('loading subjects'' times')
load /data1/intext/model8times

end

mainResDir = ['RESULTS/model' model];
disp(['Tor: making main directory' fmriDIR '/' mainResDir])
eval(['!mkdir ' fmriDIR '/' mainResDir]);
   
for JJ = subs
   snum = num2str(JJ);
   
   % make the data directory, if necessary
   % ====================================================================== 
   fmriDIR  = '/data1/intext';
	if JJ > 7,fmriDIR = '/data2/intext';,disp(['Tor: making main directory' fmriDIR '/' mainResDir]),eval(['!mkdir ' fmriDIR '/' mainResDir]);,end
	if JJ > 9,fmriDIR = '/data4/intext';,disp(['Tor: making main directory' fmriDIR '/' mainResDir]),eval(['!mkdir ' fmriDIR '/' mainResDir]);,end
	f99_CWD('LOG');                   %make the LOG and RESULTS directory 
   
if dostat
   
   % remove old files
   try 
      %eval(['cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
      disp(['!rm ' fmriDIR '/' mainResDir filesep 'sub' snum filesep '*']);
      eval(['!rm ' fmriDIR '/' mainResDir filesep 'sub' snum filesep '*']);
      disp('Removed old files')
   catch
      disp(['!rm ' fmriDIR '/' mainResDir filesep 'sub' snum filesep '*']);
   		error('Can''t find the subject result dir or erase old files')
   end
      
         % misc directory names and stuff
   % ====================================================================== 
   varname = ['sub' snum 'model8']   
   kulRES = (['RESULTS/model' model '/sub' snum]);
   statcalcRES = kulRES;
	imgdir = (['sub' snum '/task']);
   imgname = 'sra*img';
   
   % load movement parameters for user-specified covariates
   % ====================================================================== 
   eval(['cd ' fmriDIR '/' imgdir])
   P = getfiles('real*txt');
   real = load(P{1})';
   %for i = 1:6, userspec{i} = real(i,:);,end
  
      
   % make the list of directories in which to find images
   % ====================================================================== 
   ImgDir = []; ImgName = [];
   for KK = 1:numscans
      ImgDir = [ImgDir;[imgdir '/scan' num2str(KK)]];
   		ImgName = [ImgName;imgname];   
   end
   
end

   % specify and make the results directory
   % ====================================================================== 
   ResDir = ['model' model '/sub' snum];
   disp(['Tor: making data directory ' fmriDIR '/RESULTS/' ResDir]); 
   eval(['!mkdir ' fmriDIR '/RESULTS/' ResDir]);
   
   
 if dostat  
   
% calculation of all timing parameters 
% ====================================================================== 
for KK = 1:numscans

   eval(['onset{KK} = ' varname '.offset{KK};'])
	varab{KK} = zeros(size(onset{KK})); 

end


% not used - for epoch design
%length{1} = [10 0 0]; 
%length{2} = [15 0 0]; 
length = [];


% additional parameters 
% ====================================================================== 
names = ['cntylocy';'cntnlocy';'cntylocn';'cntnlocn';'probe1  ';'probe2  ']; 
TR = 2;
scanspersess = [150 150 150 150 150 150 150 150];
numconds = 6;

BasisFunct = ['ev';'ev';'ev';'ev';'ev';'ev';'ev';'ev'];
HrfType = [6 6 6 6 6 6];


% user defined hrf function 
%AA = [1:147]; 
%BB = sin(7*pi\AA); 
  

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
    TR,scanspersess,0,... 
    numconds,names,... 
    ResDir) 
  

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
    [],... 
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
    BasisFunct,... 
    [],... 
    HrfType,... 
    0,0,... 
    0,1,0,length,... 
    0,... 
    real',... 
    {'x  ';'y  ';'z  ';'rol';'pit';'yaw'}) 

% ====================================================================== 
% 'construction of the design matrix' 
% ====================================================================== 
f99_make_DesMtx 
  
disp('TorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTorTor')
disp('Construction of the design matrix finished.')
disp(['kulRES = ' kulRES])
disp(['Starting estimation'])
disp('...........................................................................')

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
    ImgDir,... 
    ImgName,... 
    statcalcRES,... 
    'None',0,... 
    'specify',[100 100 100 100 100 100 100 100],... 
    'hrf',[],... 
    'none',... 
    0)
 
 end % end if dostat
end % end loop through subjects

end % end dostat

disp('============================ tor setup ==================================')
%dostat = 0;
%docontrasts = 0;
%donorm = 0;
%dorfxcopy = 0;
%checkindss = 0;
%runfixedfx = 1;
%model = '13';
%numscans = 8;
%subs = [1 2 4 6 8 9 11 12];
%ss = subs; modelnum = model;
%drive = [1 1 1 1 2 2 4 4];


if docontrasts
   
   clear a, clear b, clear c,clear contrast, clear names

% define contrast matrices
% ==============================================================================
 a = [diag(ones(4,1),0) zeros(4,2)];
 b = zeros(4,6); 
 c{1} = [a b b b b b b b;b a b b b b b b;b b a b b b b b;b b b a b b b b;b b b b a b b b; b b b b b a b b;b b b b b b a b;b b b b b b b a];
 names{1} = 'all effects of interest         ';
 type{1} = 'F';
 
   a = [1 1 -1 -1 0 0];
  c{2} = [a a a a a a a a];
  names{2} = 'location switching              ';
  type{2} = 'T';
  
     a = [1 -1 1 -1 0 0];
  c{3} = [a a a a a a a a];
  names{3} = 'counter switching               ';
  type{3} = 'T';
  
  a = [3 -1 -1 -1 0 0 0];
  c{4} = [a a a a a a a a];
  names{4} = 'count and loc switch together   ';
  type{4} = 'T';

 c{5} = [a a a a a a a a];
 names{5} = 'fx of interest fixed across sess';
 type{5} = 'F';
 
 a = [1 1 -1 -1 0 0];
 b = [-1 -1 1 1 0 0];
 c{6} = [a a a a b b b b];
 names{6} = 'location switching x session    ';
 type{6} = 'T';

 a = [1 -1 1 -1 0 0];
 b = [-1 1 -1 1 0 0];
 c{7} = [a a a a b b b b];
 names{7} = 'counter switching x session     ';
 type{7} = 'T';
 % ==============================================================================

 
   tor_script_contrasts_modelX
end
if donorm
   normsmooth_modelX
end
if dorfxcopy
   copyresultsforrfx_modelX
end
if checkindss
	connum = '0003';
   checkindss_modelX
end
if runfixedfx
   runfixedfx_modelX
end

