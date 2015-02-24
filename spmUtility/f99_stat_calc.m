function f99_stat_calc(DIR,EXP,LULRES2,glob,BM,cLF,HParam_in,cHF,LParam_in,cVi,bFcon)

% DIR            - [matrix of directories]
% EXP            - [matrix of experiment names]
% RES            - directory where to find SPM_fMRIDesMtx
%                  directory where to write the results
% glob           - global normalisation ('Scaling' or 'None')
% BM             - burst mode? (1=yes, 0=no)
% cLF            - apply Hipgh Pass Filter? ('none' or 'specify')
% LF_nst         - specify your own cut off period? (1=yes, 0=no)
% HParam_in      - user defined cut off period for each session [1*nsess]
% cHF            - Low Pass Filter? ('none', 'Gaussian' or 'hrf')
% LParam_in      - user defined Gausssian FWHM (secs)
% cVi            - intrinsic correlations ('none' or 'AR(1)')
% bFcon          - Setup trial-specific F-contrasts? (1=yes, 0=no)


global fmriDIR fmriTEST kulRES

D = [];
for i = 1:size(DIR,1)
   d = [deblank(fmriDIR) filesep deblank(DIR(i,:))];
   D = [D; d];
end
DIR = D;

if f99_exist(fmriDIR,'RESULTS') == 0
   str = ['!mkdir ' fmriDIR  filesep 'RESULTS'];
   eval(str);
end
if f99_exist(fmriDIR,'LOG') == 0
   str = ['!mkdir ' fmriDIR filesep 'LOG'];
   eval(str);
end
if f99_exist([fmriDIR filesep 'LOG'],['STAT_' LULRES2 '.OK'])
   str = ['!rm -f ' fmriDIR filesep 'LOG' filesep 'STAT_' LULRES2 '.OK'];
   eval(str);
end


% Start logging
%--------------
logfile=[fmriDIR filesep 'LOG' filesep 'analysis.log'];
tijd=spm('time');
disp('');
disp('*****************************************************************');
f99_log(logfile,['STATISTICS  Script started at ' mat2str(tijd)]);
str = ['Start logging in ' logfile];
f99_log(logfile,str);


% checking for images
% -------------------
ns = size(DIR,1);
sess = zeros(1,ns);
pp = [];
for s=1:ns
  p = f99_P(DIR(s,:),EXP(s,:));
  nf  = size(p,1);
  str = ['  SESSION ' spm_str_manip(DIR(s,:),'t') ' : with number of timepoints = ' mat2str(nf)];
  f99_log(logfile,str);
  f99_checkP(DIR(s,:),EXP(s,:),0);
  if s==1,
    pp = p;
  else
    pp=str2mat(pp,p);
  end
  sess(s) =  size(p,1);
end
sessions{i} = cumsum(sess);
%P{i} = pp;


% checking directories
% --------------------
str = ['cd ' fmriDIR filesep 'RESULTS/' kulRES];
str
eval(str);
% load KUL_stat_gen

f99_log(logfile,' ');
f99_log(logfile,['  SPM_fMRIDesMtx.mat is loaded at   ' fmriDIR filesep ...
      kulRES]);
f99_log(logfile,['  Results are written in            ' fmriDIR filesep ...
      LULRES2]);

if ~strcmp(kulRES,LULRES2)
    msg1 = ['  *** WARNING : Results will not be written in the'];
    msg2 = ['                directory where the design matrix is loaded!!!'];
    disp(msg1);
    disp(msg2);
end


% logging analysis parameters
% ---------------------------
f99_log(logfile,' ');
f99_log(logfile,'PARAMETERS FOR ANALYSIS ');
if strcmp(lower(glob),'none')
    f99_log(logfile,'  No global normalisation')
elseif strcmp(lower(glob),'scaling')
    f99_log(logfile,'  Global normalisation is ''Scaling''');
end


if strcmp(lower(cLF),'none')
    f99_log(logfile,'  No High Pass Filter specified');
elseif strcmp(lower(cLF),'specify')
    f99_log(logfile,['  High Pass Filter is : 'mat2str(HParam_in)]);
elseif strcmp(lower(cLF),'auto')
    f99_log(logfile,'  Default High Pass Filter (based on peristimulus time)');
end

if strcmp(lower(cHF),'none')
    f99_log(logfile,'  No Low Pass Filter')
elseif strcmp(lower(cHF),'gaussian')
    f99_log(logfile,['  Low Pass Filter is ''Gaussian'' with params ' mat2str(L_Param_in)]);
elseif strcmp(lower(cHF),'hrf')
    f99_log(logfile,'  Low Pass Filter is ''hrf''');
end
  
if strcmp(lower(cVi),'none')
    f99_log(logfile,'  No Intrinsic correlations modelled');
elseif strcmp(lower(cVi),'AR(1)')
    f99_log(logfile,'  Intrinsic correlations modelled with ''AR(1)''');
 elseif strcmp(lower(cVi),'custom') 
    f99_log(logfile,'  Custom intrinsic correlation matrix (Vi)');
 end
f99_log(logfile,' ');

DIRe=DIR;

if fmriTEST==0
    kulRES = LULRES2;
    f99_spm_fmri_calc(DIRe,EXP,kulRES,glob,BM,cLF,HParam_in,cHF,LParam_in,cVi,bFcon)
    statusfile = [fmriDIR filesep 'LOG' filesep 'STAT_' LULRES2 '.OK'];
    f99_log(statusfile,'');
end








