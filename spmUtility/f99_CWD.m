function f99_CWD(DIR)
%-----------------------------------------
global fmriDIR fmriTEST;

% Check if root of working directory  exists
%-----------------------------------
if f99_exist(fmriDIR,'.') == 0
  error(['*** CRITICAL ERROR *** fmriDIR ' fmriDIR ' does not exists']);
end


% Change the working directory
%------------------------------
if nargin
  cwd = [deblank(fmriDIR) filesep deblank(DIR)];
  if f99_exist(fmriDIR,DIR) == 0
    str = ['!mkdir ' cwd];
    eval(str);
  end
else
  cwd = [deblank(fmriDIR)];
end
str = ['cd ' cwd];
eval(str);
CWD = pwd;

% Start ...logging
%--------------
disp('');
disp('*****************************************************************');
disp('* Automated SPM99 scripts - Dept. Radiology, KULeuven, Belgium  *');
disp('* http://www.kuleuven.ac.be/radiology/Research/fMRI/            *');   
disp('*****************************************************************');
logfile = [fmriDIR filesep 'LOG' filesep 'cwd.log'];
str = ['The Work Directory will be ' cwd];
f99_log(logfile,str);

disp('-----------------------------------------------------------------');
disp('');


