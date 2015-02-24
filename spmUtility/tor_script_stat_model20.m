%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Obligatory first part of the automation script                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  fmriTEST = 1 to run the simulation mode 
%           = 0 to run the actual post-processing 
%  fmriDIR  = the directory in which all your data is situated 
%----------------------------------------------------------------------- 

global fmriTEST fmriDIR;             %set global variabels 
fmriTEST = 0; 
fmriDIR  = '/data4/intext2';
f99_CWD('LOG');                   %make the LOG and RESULTS directory 

d = clock;
eval(['diary tor_script_stat_model14_' num2str(d(3)) '_' num2str(d(2)) '_' num2str(d(1))]) 

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
% model 14: evolved - my design estimation script, same as 13evoved, my Vi from each ss.
% model 15: same as before, but on ravols.
% model 18: same as 13,on intext2 newly realigned./preproc data.
%            hrf+td, mvmt params, smoothing, on ravols, n/s later.
% model 20: add 180 ms instead of 250, do on snravols - otherwise same as 18.
% model 21: use hrf alone instead of hrf + td.

disp('============================ tor setup ==================================')
dostat = 1;
docontrasts = 1;
donorm = 1;
dorfxcopy = 1;
checkindss = 0;
runfixedfx = 0;
model = '20';
numscans = 8;
totalimgs = 1200;
subs = [1 2 4 6 8 9 10 11 12];
ss = subs; modelnum = model;
drive = [4 4 4 4 4 4 4 4 4];

removeoldstat = 0;

clear TOR

% This script now includes:
%		entering movement params into model
%		custom Vi estimate
%		non-default reference slice and time bins (in spm_defaults)
%		custom regressor shift by x seconds.


% REMEMBER TO ADD HOWEVERY MANY MS TO EVENT TIMES TO CORRECT FOR SLICE TIMING SHIFT!
% programmed in below.  subtract 1 because 1st slice is 1, shift over by refslice - 1
TOR.refslice = 14;
TOR.fmriT = 28;
trsToAdd = (TOR.refslice - 1) ./ TOR.fmriT;

% ADD user-specified time, in case stim presentation time doesn't equal start of event of interest.
TOR.RT=2;
adds = 0.180;
trsToAdd = trsToAdd + (adds / TOR.RT);

mainResDir = ['RESULTS/model' model];
disp(['Tor: making main directory' fmriDIR '/' mainResDir])
eval(['!mkdir ' fmriDIR '/' mainResDir]);
   
for JJ = 1:length(subs)
   snum = num2str(subs(JJ));
   dnum = num2str(drive(JJ));
   
   % make the data directory, if necessary, when using multiple data drives
   % ====================================================================== 
   %fmriDIR  = '/data1/intext';
	%if subs(JJ) > 7,fmriDIR = '/data2/intext';,disp(['Tor: making main directory' fmriDIR '/' mainResDir]),eval(['!mkdir ' fmriDIR '/' mainResDir]);,end
	%if subs(JJ) > 9,fmriDIR = '/data4/intext';,disp(['Tor: making main directory' fmriDIR '/' mainResDir]),eval(['!mkdir ' fmriDIR '/' mainResDir]);,end
	%f99_CWD('LOG');                   %make the LOG and RESULTS directory 
   
   
   % specify and make the results directory
   % ====================================================================== 
   ResDir = ['model' model '/sub' snum];
   disp(['Tor: making data directory ' fmriDIR '/RESULTS/' ResDir]); 
   eval(['!mkdir ' fmriDIR '/RESULTS/' ResDir]);
   imgdir = (['sub' snum '/task']);

   
   if dostat
      
   % load movement parameters for user-specified covariates
   % ====================================================================== 
   try
      eval(['cd ' fmriDIR filesep imgdir]);
   catch 
      error(['COULD NOT CHANGE DIR : cd ' fmriDIR '/' imgdir])
   end
      
   P = getfiles('real*txt');
   real = load(P{1})';
   %for i = 1:6, userspec{i} = real(i,:);,end
   
   if removeoldstat
   % remove old files
   % ====================================================================== 
   disp('REMOVING OLD FILES')
   try 
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
   
	end % if removeoldstat
   
   % load the event times and format them
   % ====================================================================== 
   clear a,clear b,clear c  
	disp('loading subjects'' times')
	load /data1/intext/model8times
   varname = ['sub' snum 'model8'];
   eval(varname)
   eval(['a = ' varname '.offset;'])
   
   % add trsToAdd to event times to adjust for slice timing
   % ------------------------------------------------------
   disp(['adding ' num2str(trsToAdd) ' TRs to all event times to adjust for non-default fmriT and fmri_T0'])
   disp('--------------------------------------------------------------------------------------------------')
   disp('old:')
   a{end}(1:5,1:5)
   for whichO = 1:length(a)
      a{whichO}(a{whichO}>-.5) = a{whichO}(a{whichO}>-.5) + trsToAdd;
   end
   
   disp('new:')
   a{end}(1:5,1:5) 
   disp('Look in TOR.c to check SPM regressors entered into analysis.')
   disp('--------------------------------------------------------------------------------------------------')

   transevttimes			% makes the c variable from event times.
   
	end % if dostat



TOR	.MT=4;
% TOR.RT=2; specified above
TOR.drive = dnum; 
TOR.imgdir = [fmriDIR filesep imgdir];
TOR.resdir = [fmriDIR filesep 'RESULTS' filesep ResDir];
TOR.imgnames = 'snra*img';
TOR.subj = snum;   

if dostat
   
TOR.vardurations=0;   
TOR.Ptype='none'; 
TOR.Etype='linear';   
TOR.Rov='events';  
TOR.Cov=2;		% 1 is hrf alone, 2 with temp der.
TOR.c=c ;
TOR.v = 6;		% number of conditions/trial types in each scan
TOR.Cname = {'scsl','ncsl','scnl','ncnl','tor Probe1','tor Probe2'};	% names for each trial type
TOR.Volterra = 0;
TOR.userc = 0;
TOR.Global = 'None';
TOR.cLF = 'specify';
TOR.HParam = [100 100 100 100 100 100 100 100];
TOR.LFstr = 'none';										% high-pass filter text string - not used? - just a text message...
TOR.cHF = 'hrf';										% low-pass filter option
TOR.cVi = 'none';										% none or 'AR(1)'
TOR.bFcon = 0;
TOR.bEstNow = 1;

TOR.usregs = real';
TOR.usna = {'x  ','y  ','z  ','rol','pit','yaw'};

   
% autocorrelation matrix
% ------------------------------------
eval(['load /data2/intext/irfscans/s' snum 'avgxc avgxc'])
	% autocorr. function variable in this file for each subject is called avgxc.
% TOR.Vi = getv('make',avgxc,totalimgs);

end

% for other non-stat functions - list of directory names for each subject
TOR.allresdirs{JJ} = [fmriDIR filesep 'RESULTS' filesep ResDir];

if dostat
   
	TOR     

   
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

 a = [1 0 -1 0 1 0 -1 0 0 0 0 0];
 b = [0 1 0 -1 0 1 0 -1 0 0 0 0];
 c{7} = [a a a a a a a a;b b b b b b b b];
 names{7} = 'counter switching hrf or td     ';
 type{7} = 'F';
 
 a = [1 0 1 0 -1 0 -1 0 0 0 0 0];
 b = [-1 0 -1 0 1 0 1 0 0 0 0 0];
  c{8} = [a a a a b b b b];
  names{8} = 'location switching by half      ';
  type{8} = 'T';
  
  a = [1 0 -1 0 1 0 -1 0 0 0 0 0];
  b = [-1 0 1 0 -1 0 1 0 0 0 0 0];
  c{9} = [a a a a b b b b];
  names{9} = 'counter switching by half       ';
  type{9} = 'T';
  
  
 % ==============================================================================

 
   tor_script_contrasts_modelX
end

if donorm
   normsmooth_modelX
   
   P = str2mat(['/data3/biman/spm99/spm99/templates/T1.img']);
   P2 = P;
   
	for JJ = 1:length(subs)
   		snum = num2str(subs(JJ));
   		dnum = num2str(drive(JJ));
   
      
      datapath = [fmriDIR '/RESULTS/model' model '/sub' snum];
      % datapath = [studypath filesep 'sub' num2str(snum) filesep 'Anatomy'];
		P = strvcat(P,str2mat([datapath filesep 'sncon_0002.img']));
		P2 = strvcat(P,str2mat([datapath filesep 'sncon_0003.img']));

	end

spm_check_registration(P);
spm_orthviews('Interp',0);
spm_print
spm_check_registration(P2);
spm_orthviews('Interp',0);
spm_print

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

diary off