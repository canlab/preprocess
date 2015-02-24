function f99_header_edit(OP,DIR,EXP,rmflag)
% Smoothing or convolving
%___________________________________________________________________________
%
% Sets the origin of the images for SPM alignment 
%
% Uses:
%
% 
%
% Inputs
%
% New origin in the form of [x y z]
% *.img conforming to SPM data format (see 'Data')
%
% Outputs
%
% 
%
%__________________________________________________________________________
% @(#)spm_smooth_ui.m	2.4 99/04/01
global fmriTEST fmriDIR

% fmri init stuff
%----------------
if nargin == 3
   rmflag = '';
elseif nargin == 4
   %ok
else
   error('This script needs 3 arguments');  
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
logfile = [fmriDIR filesep 'LOG' filesep 'header_edit.log'];
tijd=spm('time');
f99_log(logfile,['Header Edit Script started at ' mat2str(tijd)]);
str = ['Start logging in ' logfile];
f99_log(logfile,str);
f99_log(logfile,'FIXING HEADERS');
f99_log(logfile,['  new origin: ' mat2str(OP)]);
f99_log(logfile,'  on image(s) ');

% get filenames and kernel width
%----------------------------------------------------------------------------
%SPMid = spm('FnBanner',mfilename,'2.4');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','HDR Edit');
spm_help('!ContextHelp','spm_header_edit.m');

%s     = spm_input('smoothing {FWHM in mm}',1);
origin = ['ORIGIN =' mat2str(OP) ';'];

% cycle over sessions
%--------------------
for i = 1:size(EXP)
   % removal of  statusfile
   % -----------------------------------------
   if f99_exist([fmriDIR filesep 'LOG'],'HEADER.OK')
      str = ['!rm -f ' fmriDIR filesep 'LOG' filesep 'HEADER.OK'];
      eval(str);
   end
   
   % setup and control images
   % ---------------------------------------------------------------------------
   f99_log(logfile,['   ' mat2str(DIR(i,:)) filesep mat2str(EXP(i,:))]);
   f99_checkP(DIR(i,:),EXP(i,:),0);
   
   str = [DIR(i,:) filesep EXP(i,:)];
   [path,file,ext] = fileparts(str);
   str = [path filesep file '.mat'];
   matExist = dir(str);

   if size(matExist,1) ~= 0
      f99_log(logfile,'     WARNING ! mat Files Already Exist');
      if isempty(rmflag)
      	f99_log(logfile,'     WARNING ! Header Edit will not take effect due to existence of mat files!');
      else
      	f99_log(logfile,'     WARNING ! existing mat files will be removed!');
      end
   end
   
   %P     = spm_get(Inf,'.img','select scans');
   P = f99_P(DIR(i,:),EXP(i,:));
   n     = size(P,1);
   
   % implement the header change
   %---------------------------------------------------------------------------
   if fmriTEST == 0
   
    statusfile = [fmriDIR filesep 'LOG' filesep 'HEADER.OK'];
    f99_log(statusfile,'');
    if ~isempty(rmflag)
    	str = [DIR(i,:) filesep EXP(i,:)];
   	[path,file,ext] = fileparts(str);
   	str = [file '.mat'];
    	f99_rm_file(rmflag,DIR(i,:),str);
    end
    
      spm('Pointer','Watch');
      spm('FigName','Header Edit: working',Finter,CmdLine);
      spm_progress_bar('Init',n,'Header Edit','Volumes Complete');
      for j = 1:n
         Q = deblank(P(j,:));
         [pth,nm,xt,vr] = fileparts(deblank(Q));
         U = fullfile(pth,[nm '.hdr' vr]);
         status = spm_header_edit('set',U,origin);
         if status == 1
         	f99_log(logfile,['     WARNING ! header edit failure for file:' U]);
         end
         
         spm_progress_bar('Set',j);
      end
      spm_progress_bar('Clear',j);
      spm('FigName','Header Edit: done',Finter,CmdLine);
      spm('Pointer');
    end  
    
  end
  