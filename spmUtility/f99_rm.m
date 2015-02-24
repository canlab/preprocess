function f99_rm(check,DIR,EXP,short)
%==============================================
% check - als deze status file bestaat
% DIR - matrix van directories
% EXP - matrix van experimenten
%  van de files die moeten verwijderd worden
%
% optioneel:
% short 
%   0 = lange namen (met dir ervoor)
%   1 = korte namen, met gebruik van fmriDIR 
%----------------------------------------------
global fmriTEST fmriDIR;

% fmri init stuff
%----------------
if nargin == 3
  short = 0;
end

% fmri logging stuff
%-------------------
logfile = [fmriDIR filesep 'LOG' filesep 'remove.log'];

if short == 1
  D = [];
  for i = 1:size(DIR,1)
    d = [deblank(fmriDIR) filesep deblank(DIR(i,:))];
    D = [D; d];
  end
  DIR = D;
end

str = [fmriDIR filesep 'LOG'];
if f99_exist(str,check) == 1
  for i = 1:size(DIR,1)
    if findstr('*',DIR(i,:))
      f99_log(logfile,'     WARNING no wildcards accepted in directories !');
    else
      if strcmp('*',EXP(i,size(EXP(i,:),2))) == 1 | ~isempty(findstr(' *',EXP(i,:))) | ~isempty(findstr('* ',EXP(i,:)))
        f99_log(logfile,'     WARNING no wildcards with _*_ accepted in experiments !');
      elseif ~isempty(findstr('anat',DIR(i,:)))
        f99_log(logfile,'     WARNING no remove of anatomy files !');
      else
        if fmriTEST
          f99_log(logfile,['     will rm ' deblank(DIR(i,:)) filesep deblank(EXP(i,:)) ]);
        else
          str = ['!rm -f ' deblank(DIR(i,:)) filesep deblank(EXP(i,:)) ];
          str2 = spm_str_manip(str,'s');
          str2 = [deblank(str2) '.hdr'];
          f99_log(logfile,str);
          f99_log(logfile,str2);
          eval(str)
          eval(str2)
        end
      end
    end
  end
else
%  f99_log(logfile,['     if ' check ' exists then']);
  for i = 1:size(DIR,1)
    if findstr('*',DIR(i,:))
      f99_log(logfile,'     WARNING no wildcards accepted in directories !');
    else
      if strcmp('*',EXP(i,size(EXP(i,:),2))) == 1 | ~isempty(findstr(' *',EXP(i,:))) | ~isempty(findstr('* ',EXP(i,:)))
        f99_log(logfile,'     WARNING no wildcards with _*_ accepted in experiments !');
      elseif ~isempty(findstr('anat',DIR(i,:)))
        f99_log(logfile,'     WARNING no remove of anatomy files !');
      else
        f99_log(logfile,['     if ' check ' exists then will REMOVE these files']);
      end
    end
  end
end

