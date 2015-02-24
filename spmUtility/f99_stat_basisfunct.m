function f99_stat_basisfunct(Rov,hrffunct,Cov,pst_ev,h,m,HRF,TD,W,VOLT,regr_us,regr_us_na)

% 'basisfunctions to use'
% Rov         - matrix with types of trials 
%                * event ('ev') or epoch ('ep')
%                * one string per trial type (e.g. ['ep';'ep';'ev'])
% hrffunct    - user specified hrf function of this subject
%                length of this function can be checked using fmriTEST
% Cov         - model type for covariates of interest
%                1 = basis functions  (Discrete Cosine Set)
%                2 = basis functions  (Mean & exponential decay)
%                3 = fixed response   (Half-sine)
%                4 = fixed response   (Box-car)
%             - model type for event-related responses
%                5 = hrf (alone)
%                6 = hrf (with time derivative)
%                7 = hrf (with time and dispersion derivatives)
%                8 = basis functions (Fourier set)
%                9 = basis functions (Windowed Fourier set)
%               10 = basis functions (Gamma functions)
%               11 = basis functions (Gamma functions with
%                    derivatives)
%                * e.g. [4 4 6]
% 'event'
% pst_ev      - window length          (for Cov=8 Cov=9)
% h           - order of Fourier set   (for Cov=8 Cov=9)
% 'epoch'
% m           - number of basis functions (if Cov=1)
% HRF         - convolve with HRF (1=yes / 0=no)
% TD          - add temporal differences (1=yes / 0=no)
% W           - epoch length (scans) for each condition  
% 'additional'
% VOLT        - interactions among trials (Volterra series)? (1=yes / 0=no)
% regr_us     - matrix with user specified regressors
% regr_us_na  - matrix with strings of names of regr_us

global fmriTEST fmriDIR kulRES
logfile=[fmriDIR filesep 'LOG' filesep 'design.log'];

f99_log(logfile,['BASISFUNCTIONS']);

str_1  = 'basis functions (Discrete Cosine Set)';
str_2  = 'basis functions (Mean & exponential decay)';
str_3  = 'fixed response (Half-sine)';
str_4  = 'fixed response (Box-car)';
str_5  = 'hrf (alone)';
str_6  = 'hrf (with time derivative)';
str_7  = 'hrf (with time and dispersion derivatives)';
str_8  = 'basis functions (Fourier set)';
str_9  = 'basis functions (Windowed Fourier set)';
str_10 = 'basis functions (Gamma functions)';
str_11 = 'basis functions (Gamma functions with derivatives)';


for i = 1:length(Cov)
  if strcmp(deblank(Rov(i,:)),'ev')
    st_Rov = 'event';
  elseif strcmp(deblank(Rov(i,:)),'ep')
    st_Rov = 'epoch';
  else
    error('You typed in a wrong Trial Type (should be: ''ev'' or ''ep'')');
  end
  str1 = ['  Trial ' mat2str(i) ' is ' st_Rov ','];
  strP = ['str = str_' mat2str(Cov(i)) ' ;'];
  eval(strP);
  str2 = [' modelled with ' str];
  f99_log(logfile,[str1 str2]);
end

filename = [fmriDIR filesep 'RESULTS' filesep kulRES filesep 'KUL_stat_basisfunct'];
str = ['save ' filename ' Rov hrffunct Cov pst_ev h m HRF TD W VOLT regr_us regr_us_na'];
eval(str)

