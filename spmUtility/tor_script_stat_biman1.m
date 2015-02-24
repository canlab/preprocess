%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Obligatory first part of the automation script                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  fmriTEST = 1 to run the simulation mode 
%           = 0 to run the actual post-processing 
%  fmriDIR  = the directory in which all your data is situated 
%----------------------------------------------------------------------- 

global fmriTEST fmriDIR;             %set global variabels 
fmriTEST = 0; 
fmriDIR  = '/data/biman4';
f99_CWD('LOG');                   %make the LOG and RESULTS directory 

d = clock;
eval(['diary tor_script_stat_biman1_' num2str(d(3)) '_' num2str(d(2)) '_' num2str(d(1))]) 

% THIS SCRIPT NOW COREGISTERS AND SEGMENTS T1s TOO (if normalizing)! CAREFUL!

% model1: temporal smoothing, no HP filter


disp('============================ tor setup ==================================')
dostat = 0;
removeoldstat = 0;
docontrasts = 1;
donorm = 0;
dorfxcopy = 0;
checkindss = 0;
runfixedfx = 0;
model = '1';
numscans = 6;
totalimgs = 1632;
subs = [3 4 5 6 8 9 10];
ss = subs; modelnum = model;
drive = [''];


clear TOR
TOR.numscans = numscans;
TOR.evttimesfile = [fmriDIR filesep 'bimanevttimes.mat'];	% name of file with event times
													% format: sub*.offset{*} is scan
													% in this cell, rows are conditions, cols are list of times
													% padded with -1's.

% This script now includes:
%		entering movement params into model
%		custom Vi estimate
%		non-default reference slice and time bins (in spm_defaults)
%		custom regressor shift by x seconds.


% REMEMBER TO ADD HOWEVERY MANY MS TO EVENT TIMES TO CORRECT FOR SLICE TIMING SHIFT!
% programmed in below.  subtract 1 because 1st slice is 1, shift over by refslice - 1
TOR.refslice = 17;
TOR.fmriT = 35;
trsToAdd = (TOR.refslice - 1) ./ TOR.fmriT;

global fMRI_T
global fMRI_T0

fMRI_T = TOR.fmriT;
fMRI_T0 = TOR.refslice;

% ADD user-specified time, in case stim presentation time doesn't equal start of event of interest.
TOR.RT=2.5;
adds = 0;	% user-specified regressor shift from stimulus onset.
trsToAdd = trsToAdd + (adds / TOR.RT);

mainResDir = ['RESULTS/model' model];

% eval(['if ~(exist(''' fmriDIR filesep mainResDir ''') == 7),!mkdir ' fmriDIR '/' mainResDir ', end']);

myMainResDir = [fmriDIR filesep mainResDir];
if ~(exist(myMainResDir)) == 7
	disp(['Tor: making main directory ' myMainResDir])
	eval(['!mkdir ' myMainResDir])
end




if donorm
   spm fmri
   coregnormpath
   % normsmooth_modelX done in coregnormpath

end












for JJ = 1:length(subs)
   snum = num2str(subs(JJ));
   if ~isempty(drive),dnum = num2str(drive(JJ));,else dnum = [];,end
   
   % make the data directory, if necessary, when using multiple data drives
   % ====================================================================== 
   %fmriDIR  = '/data1/intext';
	%if subs(JJ) > 7,fmriDIR = '/data2/intext';,disp(['Tor: making main directory' fmriDIR '/' mainResDir]),eval(['!mkdir ' fmriDIR '/' mainResDir]);,end
	%if subs(JJ) > 9,fmriDIR = '/data4/intext';,disp(['Tor: making main directory' fmriDIR '/' mainResDir]),eval(['!mkdir ' fmriDIR '/' mainResDir]);,end
	%f99_CWD('LOG');                   %make the LOG and RESULTS directory 
   
   
   % specify and make the results directory
   % ====================================================================== 
   ResDir = ['model' model '/sub' snum];

   myResDir = [fmriDIR '/RESULTS/' ResDir];
   if ~(exist(myResDir)) == 7
   	disp(['Tor: making data directory ' myResDir]); 
	eval(['!mkdir ' myResDir])
   end

   %eval(['if ~(exist(''' fmriDIR '/RESULTS/' ResDir ''') == 7),!mkdir ' fmriDIR '/RESULTS/' ResDir ', end']);
   
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
   	end % if removeoldstat
   
   	try
		eval(['cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
	catch
      		error(['CANNOT CHANGE DIR TO : cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
   	end

   	% load the event times and format them
   	% ====================================================================== 
   	clear a,clear b,clear c  
	
	str = (['load ' TOR.evttimesfile]);
	disp(['loading subjects'' times with command: ' str])
	eval(str)

   	varname = ['sub' snum];
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



   TOR.MT=4;
   TOR.drive = dnum; 

   
   for k = 1:numscans
	TOR.imgdir{k} = [fmriDIR filesep imgdir filesep 'scan' num2str(k)];
   end

   TOR.resdir = [fmriDIR filesep 'RESULTS' filesep ResDir];
   TOR.imgwildcard = 'snra*img';

   % Build list of image names
   TOR.imgnames = tor_list_files(TOR.imgdir,TOR.imgwildcard);

   TOR.subj = snum;   

   if dostat
   
	TOR.vardurations=0;   
	TOR.Ptype='none'; 
	TOR.Etype='linear';   
	TOR.Rov='events';  
	TOR.Cov=1;		% 1 = hrf alone, 2 = hrf + td
	TOR.c=c ;
	TOR.v = 10;		% number of conditions/trial types in each scan
	TOR.Cname = {'IL','IR','CL','CR','SR','SL','CR','CL','errors','instructions'};	% names for each trial type
	TOR.Volterra = 0;
	TOR.userc = 0;
	TOR.Global = 'None';
	TOR.cLF = 'none';
	TOR.HParam = [0 0 0 0 0 0];								% hp filter cutoff values.
	TOR.LFstr = 'none';										% high-pass filter text string - not used? - just a text message...
	TOR.cHF = 'hrf';										% low-pass filter option - hrf, gaussian, none
	TOR.cVi = 'none';										% none or 'AR(1)'
	TOR.bFcon = 0;
	TOR.bEstNow = 1;

	TOR.usregs = real';
	TOR.usna = {'x  ','y  ','z  ','rol','pit','yaw'};

	TOR.autocorr = '';									% autocorrelation model: name of acf to load or empty for spm's way.
	
	% autocorrelation matrix
	% ------------------------------------
	if ~isempty(TOR.autocorr)
		eval(['load /data2/intext/irfscans/s' snum 'avgxc avgxc'])
		% autocorr. function variable in this file for each subject is called avgxc.
		TOR.Vi = getv('make',avgxc,totalimgs);
	end

   end	% if dostat

   % for other non-stat functions - list of directory names for each subject
   TOR.allresdirs{JJ} = [fmriDIR filesep 'RESULTS' filesep ResDir];





   if dostat
   
	TOR     

	% change to results directory!
	% ------------------------------------
	try
      		eval(['cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
	catch
      		error(['CANNOT CHANGE DIR TO : cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
	end
 


   	% build and estimate the model
   	% ====================================================================== 
   	tor_spm_fmri_spm_ui(TOR)
   
   
   end		% dostat

end		% loop thru subjects


if docontrasts
   clear a, clear b, clear c,clear contrast, clear names



for JJ = 1:length(subs)
   	snum = num2str(subs(JJ));

	
	% change to results directory!
	% ------------------------------------
	try
      		eval(['cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
		load SPM.mat
		
	catch
      		error(['CANNOT CHANGE DIR TO : cd ' fmriDIR '/' mainResDir filesep 'sub' snum]);
	end

	% just remove the old xCon.mat
		!rm xCon.mat
		!rm con*img
		!rm con*hdr
		!rm spmT*img
		!rm spmT*hdr
		!rm spmF*img
		!rm spmF*hdr

	disp('*****************************************************************')   
	disp('Defining contrasts')

% define contrast matrices
% ==============================================================================
emptyCols = sum(xX.X) == 0;

% CONTRAST 1
cnum = 1;
a = [1 1 1 1 1 1 1 1 0 0];
c{cnum} = tor_make_F_contrast_vector(a,Sess);
 names{1} = 'all effects of interest         ';
 type{1} = 'F';
 

% CONTRAST 2
   cnum = 2;
   a = [1 1 1 1 -1 -1 -1 -1 0 0];
   c{cnum} = tor_make_T_contrast_vector(a,Sess);
   names{cnum} = 'dual vs. single attribute       ';
   type{cnum} = 'T';
  
% CONTRAST 3 
cnum = 3;
  a = [1 1 -1 -1 0 0 0 0 0 0];
   c{cnum} = tor_make_T_contrast_vector(a,Sess);    
   names{3} = 'incongruent vs. congruent       ';
   type{3} = 'T';
 
% CONTRAST 4 
cnum = 4;
  a = [1 -1 1 -1 -1 1 -1 1 0 0];
  c{cnum} = tor_make_T_contrast_vector(a,Sess);
  names{4} = 'left vs. right                   ';
  type{4} = 'T';

% CONTRAST 5 
cnum = 5;
  a = [-1 1 -1 1 1 -1 1 -1 0 0];
  c{cnum} = tor_make_T_contrast_vector(a,Sess);
  names{5} = 'right vs. left                   ';
  type{5} = 'T';

% CONTRAST 6 
cnum = 6;
  a = [0 0 0 0 1 1 -1 -1 0 0];
  c{cnum} = tor_make_T_contrast_vector(a,Sess);
  names{6} = 'shape vs. color                  ';
  type{6} = 'T';
 
% CONTRAST 7 
cnum = 7;
  a = [0 0 0 0 -1 -1 1 1 0 0];
  c{cnum} = tor_make_T_contrast_vector(a,Sess);
  names{7} = 'color vs. shape                  ';
  type{7} = 'T';


% CONTRAST 8
cnum = 8;
   a = [0 0 0 0 0 0 0 0 1 0]; 
   c{cnum} = tor_make_F_contrast_vector(a,Sess);
   names{cnum} = 'errors                          ';
   type{cnum} = 'F';  
 



% CONTRAST 8
%cnum = 8;
%   a = [1 1 1 1 -1 -1 -1 -1 0 0]; b = -a;
%   a = a ./ sum(a(a>0)); b = b ./ sum(b(b>0))
%   if ~(sum(a) == 0),error('ERROR IN CONTRAST: DOES NOT SUM TO 0'),end
%   c{8} = [a a a b b b];
%   names{8} = 'dual vs. single by practice     ';
%   type{8} = 'T';
  
% CONTRAST 9
%cnum = 9;
%   a = [1 1 -1 -1 0 0 0 0 0 0]; b = -a;
%   a = a ./ sum(a(a>0)); b = b ./ sum(b(b>0))
%   if ~(sum(a) == 0),error('ERROR IN CONTRAST: DOES NOT SUM TO 0'),end
%   c{9} = [a a a b b b];
%   names{9} = 'incong vs. cong by practice     ';
%   type{9} = 'T';

% CONTRAST 11
%cnum = 11;
% c{cnum} = tor_make_F_contrast_vector(a,Sess);
% names{cnum} = 'fx of interest fixed across sess';
% type{cnum} = 'F';
 
 % ==============================================================================

 
   estimateContrasts

end % loop thru subjects

end % if docontrasts



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