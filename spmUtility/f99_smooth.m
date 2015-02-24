function f99_smooth(SP,DIR,EXP,rmflag)

% '==========================================================================='
%  ' SMOOTHING OF FUNCTIONAL IMAGES'
% '==========================================================================='
% '               --- The Script Parameters Explained ---                     '
%
%
%---------------------------------------------------------------------------
%  SP  - FWHM of smoothing kernel
%  DIR  - [matrix of dirs]
%  EXP  - [matrix of experiments]
%---OPTIONALY----------------------------------------------------------------
%   rmflags - status file in the LOG that controls if files can be
%             safely deleted
%---------------------------------------------------------------------------
%
%
%function spm_smooth_ui
% Smoothing or convolving
%___________________________________________________________________________
%
% Convolves image files with an isotropic (in real space) Gaussian kernel 
% of a specified width.
%
% Uses:
%
% As a preprocessing step to suppress noise and effects due to residual 
% differences in functional and gyral anatomy during inter-subject 
% averaging.
%
% Inputs
%
% *.img conforming to SPM data format (see 'Data')
%
% Outputs
%
% The smoothed images are written to the same subdirectories as the 
% original *.img and are prefixed with a 's' (i.e. s*.img)
%
%__________________________________________________________________________
% @(#)spm_smooth_ui.m	2.7	99/09/20

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'smooth',index_of_Analysis};
%_______________________________________________________________________

global BCH fmriTEST fmriDIR;

% fmri init stuff
%----------------
if nargin == 3
   rmflag = '';
elseif nargin == 4
   %ok
else
   error('This script needs 3 arguments (+ 1 optional rmflag)');   
end

D = [];
for i = 1:size(DIR,1)
   d = [deblank(fmriDIR) filesep deblank(DIR(i,:))];
   D = [D; d];
end
DIR = D;

if f99_exist(fmriDIR,'RESULTS') == 0
   str = ['!mkdir ' fmriDIR filesep 'RESULTS'];
   eval(str);
end
if f99_exist(fmriDIR,'LOG') == 0
   str = ['!mkdir ' fmriDIR filesep 'LOG'];
   eval(str);
end

% Start logging
%--------------
disp('');
disp('*****************************************************************');
logfile = [fmriDIR filesep 'LOG' filesep 'smoothing.log'];
tijd=spm('time');
f99_log(logfile,['Smoothing Script started at ' mat2str(tijd)]);
str = ['Start logging in ' logfile];
f99_log(logfile,str);
f99_log(logfile,'SMOOTHING');
f99_log(logfile,['  kernel ' mat2str(SP)]);
f99_log(logfile,'  on image(s) ');

% get filenames and kernel width
%----------------------------------------------------------------------------
%SPMid = spm('FnBanner',mfilename,'2.4');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Smooth');
spm_help('!ContextHelp','spm_smooth_ui.m');

%s     = spm_input('smoothing {FWHM in mm}',1,...
%                  'batch',{},'FWHMmm');
s = SP;

%if isempty(BCH)
%   P = spm_get(Inf,'.img','select scans');
%else
%   P = spm_input('batch',{},'files');
%end

% cycle over sessions
%--------------------
for i = 1:size(EXP)
   % removal of  statusfile
   % -----------------------------------------
   if f99_exist([fmriDIR filesep 'LOG'],'SMOOTH.OK')
      str = ['!rm -f ' fmriDIR filesep 'LOG' filesep 'SMOOTH.OK'];
      eval(str);
   end
   
   % setup and control images
   % ---------------------------------------------------------------------------
   f99_log(logfile,['   ' mat2str(DIR(i,:)) filesep mat2str(EXP(i,:))]);
   f99_checkP(DIR(i,:),EXP(i,:),0);
   if EXP(i,1) ~= 'n'
      f99_log(logfile,'     WARNING ! are the images normalised ????');
   end
   
   %P     = spm_get(Inf,'.img','select scans');
   P = f99_P(DIR(i,:),EXP(i,:));
   n     = size(P,1);

   % implement the convolution
   %---------------------------------------------------------------------------
   if fmriTEST == 0
     spm('Pointer','Watch');
     spm('FigName','Smooth: working',Finter,CmdLine);
     spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
     for j = 1:n
	Q = deblank(P(j,:));
	[pth,nm,xt,vr] = fileparts(deblank(Q));
	U = fullfile(pth,['s' nm xt vr]);
	spm_smooth(Q,U,s);
	spm_progress_bar('Set',j);
     end
     spm_progress_bar('Clear',j);
     spm('FigName','Smooth: done',Finter,CmdLine);
     spm('Pointer');
     end
end

statusfile = [fmriDIR filesep 'LOG' filesep 'SMOOTH.OK'];
f99_log(statusfile,'');
if ~isempty(rmflag)
        f99_rm(rmflag,DIR,EXP);
end
