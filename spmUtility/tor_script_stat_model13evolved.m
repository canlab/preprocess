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
%f99_CWD('LOG');                   %make the LOG and RESULTS directory 

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
docontrasts = 0;
donorm = 0;
dorfxcopy = 1;
checkindss = 0;
runfixedfx = 0;
model = '13';
numscans = 8;
subs = [1 2 4 6 8 9 10 11 12];
ss = subs; modelnum = model;
drive = [1 1 1 1 2 2 4 4 4];



mainResDir = ['RESULTS/model' model];
disp(['Tor: making main directory' fmriDIR '/' mainResDir])
eval(['!mkdir ' fmriDIR '/' mainResDir]);
   
for JJ = 1:length(subs)
   snum = num2str(subs(JJ));
   dnum = num2str(drive(JJ));
   
   % make the data directory, if necessary
   % ====================================================================== 
   fmriDIR  = '/data1/intext';
	if subs(JJ) > 7,fmriDIR = '/data2/intext';,disp(['Tor: making main directory' fmriDIR '/' mainResDir]),eval(['!mkdir ' fmriDIR '/' mainResDir]);,end
	if subs(JJ) > 9,fmriDIR = '/data4/intext';,disp(['Tor: making main directory' fmriDIR '/' mainResDir]),eval(['!mkdir ' fmriDIR '/' mainResDir]);,end
	%f99_CWD('LOG');                   %make the LOG and RESULTS directory 
   
   
   % specify and make the results directory
   % ====================================================================== 
   ResDir = ['model' model '/sub' snum];
   disp(['Tor: making data directory ' fmriDIR '/RESULTS/' ResDir]); 
   eval(['!mkdir ' fmriDIR '/RESULTS/' ResDir]);
   
   
if dostat
   % load movement parameters for user-specified covariates
   % ====================================================================== 
   imgdir = (['sub' snum '/task']);
   try
      eval(['cd ' fmriDIR '/' imgdir]);
   catch 
      error(['COULD NOT CHANGE DIR : cd ' fmriDIR '/' imgdir])
   end
      
   P = getfiles('real*txt');
   real = load(P{1})';
   %for i = 1:6, userspec{i} = real(i,:);,end

   
   % remove old files
   % ====================================================================== 
   disp('REMOVING OLD FILES')
   try 
      %eval(['cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
      disp(['!rm ' fmriDIR '/' mainResDir filesep 'sub' snum filesep '*']);
      eval(['!rm ' fmriDIR '/' mainResDir filesep 'sub' snum filesep '*']);
      disp('Removed old files')
   catch
      disp(['!rm ' fmriDIR '/' mainResDir filesep 'sub' snum filesep '*']);
   		error('Can''t find the subject result dir or erase old files')
   end
   try
      eval(['cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
   catch
      error(['CANNOT CHANGE DIR TO : cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
	end
   

   
   % load the event times and format them
   % ====================================================================== 
   clear a,clear b,clear c  
	disp('loading subjects'' times')
	load /data1/intext/model8times
   varname = ['sub' snum 'model8'];
   eval(varname)
   eval(['a = ' varname '.offset;'])
   transevttimes			% makes the c variable from event times.
   
   clear TOR

TOR	.MT=4;
TOR.RT=2;
TOR.drive = dnum;   
TOR.subj = snum;   
TOR.vardurations=0;   
TOR.Ptype='none'; 
TOR.Etype='linear';   
TOR.Rov='events';  
TOR.Cov=2;
TOR.c=c ;
TOR.Volterra = 0;
TOR.userc = 0;
TOR.Global = 'None';
TOR.cLF = 'specify';
TOR.HParam = [100 100 100 100 100 100 100 100];
TOR.LFstr = 'hrf';
TOR.cHF = 'hrf';
TOR.cVi = 'none';
TOR.bFcon = 0;
TOR.bEstNow = 1;

TOR.usregs = real';
TOR.usna = {'x  ','y  ','z  ','rol','pit','yaw'};


      

   
   % build and estimate the model
   % ====================================================================== 
   tor_spm_fmri_spm_ui(TOR)
   
   
end		% dostat
end		% loop thru subjects

   
   
   

if docontrasts
   
   clear a, clear b, clear c,clear contrast, clear names

% define contrast matrices
% ==============================================================================
 a = [diag(ones(8,1),0) zeros(8,4)];
 b = zeros(8,12); 
 c{1} = [a b b b b b b b;b a b b b b b b;b b a b b b b b;b b b a b b b b;b b b b a b b b; b b b b b a b b;b b b b b b a b;b b b b b b b a];
 names{1} = 'all effects of interest         ';
 type{1} = 'F';
 
   a = [1 0 1 0 -1 0 -1 0 0 0 0 0];
  c{2} = [a a a a a a a a];
  names{2} = 'location switching              ';
  type{2} = 'T';
  
     a = [1 0 -1 0 1 0 -1 0 0 0 0 0];
  c{3} = [a a a a a a a a];
  names{3} = 'counter switching               ';
  type{3} = 'T';
  
  a = [3 0 -1 0 -1 0 -1 0 0 0 0 0];
  c{4} = [a a a a a a a a];
  names{4} = 'count and loc switch together   ';
  type{4} = 'T';

 c{5} = [a a a a a a a a];
 names{5} = 'fx of interest fixed across sess';
 type{5} = 'F';
 
 a = [1 0 1 0 -1 0 -1 0 0 0 0 0];
 b = [0 1 0 1 0 -1 0 -1 0 0 0 0];
 c{6} = [a a a a a a a a;b b b b b b b b];
 names{6} = 'location switching hrf or td    ';
 type{6} = 'F';

 a = [1 0 1 0 -1 0 -1 0 0 0 0 0];
 b = [0 1 0 -1 0 1 0 -1 0 0 0 0];
 c{7} = [a a a a a a a a;b b b b b b b b];
 names{7} = 'counter switching hrf or td     ';
 type{7} = 'F';
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

