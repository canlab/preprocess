function f99_stat_modul(Ptype,Pstr,Etype,h_exp,h_pol,sel,p_mod)

% 'respons modulation'
% Ptype       - parametric modulation? ('none','time','other')
% Pstr        - name of parameter 
%                    (in case of 'other')
% Etype       - expansion type ('linear','exponen','polynom')
% h_exp       - time constant 
%                    (in case of 'exponen')
% h_pol       - order of polynomial expansion 
%                    (in case of 'polynom')
% sel         - which trials to modulate
% p           - matrix with parameters for each selected covariate
%               [#starting points per covariate * #selected covariates]

global fmriTEST fmriDIR kulRES
logfile=[fmriDIR filesep 'LOG' filesep 'design.log'];

if strcmp(deblank(lower(Ptype)),'none')
  f99_log(logfile,['NO MODULATIONS']);
  Pstr    = '';
  Etype   = '';
  h_exp   = 0;
  h_pol   = 0;
  sel     = [];
  p_mod   = [];
else
  f99_log(logfile,['MODULATIONS']);
  if strcmp(deblank(lower(Ptype)),'time')
    f99_log(logfile,['  parametric modulation with time']);
  else
    f99_log(logfile,['  parametric modulation with ' Pstr]);
  end
  if strcmp(deblank(lower(Etype)),'linear')
    f99_log(logfile,['  linear expansion']);
  end
  if strcmp(deblank(lower(Etype)),'exponen')
    f99_log(logfile,['  exponential expansion with time constant ' mat2str(h_exp)]);
  end
  if strcmp(deblank(lower(Etype)),'polynom')
    f99_log(logfile,['  polynomial expansion of order ' mat2str(h_pol)]);
  end
  f99_log(logfile,['  trials ' mat2str(sel) ' are modulated']);
  if strcmp(deblank(lower(Ptype)),'other')
    for i = 1 : size(p_mod,1)
      f99_log(logfile,['  modulation parameters for ' Pstr ' are']);
      f99_log(logfile,['       ' mat2str(p_mod(i,:))]);
    end
  end
end

filename = [fmriDIR filesep 'RESULTS' filesep kulRES filesep 'KUL_stat_modul'];
str = ['save ' filename ' Ptype Pstr Etype h_exp h_pol sel p_mod'];
eval(str)

