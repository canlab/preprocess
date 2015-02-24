function f99_stat_stoch(s_STOC,nulev,soa_stoc,oc_prob,prob_type)

% 'stochastic design'
% DO YOU WANT A STOCHASTIC DESIGN ('yes'/'no')
%    if you answer yes here, answer no in the next part
% nulev       - include null event? (1=yes / 0=no ; for each session)
% soa_stoc    - soa for each trial type
% oc_prob     - matrix with occurence probabilities for each trial type
% prob_type   - stationary (1) or modulated (0) occurence probabilities

global fmriTEST fmriDIR kulRES
logfile=[fmriDIR filesep 'LOG' filesep 'design.log'];


STOC = 1;
if strcmp(deblank(lower(s_STOC)),'no')
  nulev     = 0;
  soa_stoc  = [];
  oc_prob   = [];
  prob_type = 0;
  STOC      = 0;
else
  f99_log(logfile,['STOCHASTIC']);
  if any(nulev)
    f = find(nulev);
    f99_log(logfile,['  Null event is included in session(s) ' mat2str(f)]);
  else
    f99_log(logfile,['  NO null event is included']);
  end
  f99_log(logfile,['  Stimulus Onset Asynchronies for trial types = ' mat2str(soa_stoc)]);
  f99_log(logfile,['  Matrix wih occurence probabilities = ' mat2str(oc_prob)]);
  if any(prob_type)
    g = find(prob_type);
    f99_log(logfile,['  Probabilities are stationary in session(s) ' mat2str(g)]);
  end
  if any(prob_type==0)
    h = find(~prob_type);
    f99_log(logfile,['  Probabilities are modulated  in session(s) ' mat2str(h)]);
  end
end

filename = [fmriDIR filesep 'RESULTS' filesep kulRES filesep 'KUL_stat_stoch'];
str = ['save ' filename ' STOC nulev soa_stoc oc_prob prob_type'];
eval(str)

