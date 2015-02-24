% * BEGIN USER INPUT

% don't need these for revised TOR spm stats
% global fmriTEST fmriDIR;             %set global variabels 
% fmriTEST = 0; 
% fmriDIR  = '/data/biman5';
% f99_CWD('LOG');                   %make the LOG and RESULTS directory 

% main experiment directory
% ====================================================================== 
fmriDIR  = EXPT.SUBJECT.studydir;


if ~(exist('SubjCode') == 1)
	SubjCode = input('Enter subject code ','s');
end

% open a diary file to save output
% ====================================================================== 
diary Intext_Model1_Script_Diary


% define the vectors of onsets for each trial type
% ====================================================================== 

% For Bimanual interference task
%	1:11 IncongR IncongL CongR CongL ShapeR ShapeL ColorR ColorL Baseline Shape_Instruct Color_Instruct
%	12:14	int_blk	shape_blk col_blk	
%	c variable contains cell arrays for these event types.

cd C:\Tor_Documents\CurrentExperiments\intext2
load intext_c
eval(['c = ' SubjCode '.c;'])


% *********************************************************************************************************
% *                                                                                                       *
% *  START USER INPUT                                                                                     *
% *                                                                                                       *
% *********************************************************************************************************

modelname = 'model1';

% define input parameters to spm_spm_ui
% ====================================================================== 

TOR.estimatemodel = 1;				% estimate the statistical model
TOR.estimatecontrasts = 1;			% estimate all specified contrasts
TOR.getresults = 1;					% get and print the results to a ps file

TOR.refslice = 1;					% reference slice for slice timing.  default 1
TOR.fmriT = 16;						% time bins per TR in delta functions.  default 16
TOR.RT=2;						    % TR of your study

TOR	.MT=4;												
	% 1 = specify
	% 2 = review
	% 3 = estimate
	% 4 = specify and estimate

for i = 1:EXPT.nsess
    TOR.imgdir{i} = [fmriDIR filesep SubjCode filesep 'scan' num2str(i)];
end

TOR.resdir = [fmriDIR filesep 'RESULTS' filesep modelname filesep SubjCode];

TOR.imgwildcard = 'snra*img';

TOR.numscans = EXPT.nsess;
	%number of sessions in design	
	
TOR.vardurations=0;   
TOR.Ptype='none'; 
TOR.Etype='linear';   
TOR.Rov='events';

TOR.Cov=1;
	% For Events:	1 = hrf alone, 2 = hrf + td
	% For Epochs:	1 = 'basis functions  (Discrete Cosine Set)',...
	%			2 = 'basis functions  (Mean & exponential decay)',...
	%			3 = 	 'fixed response   (Half-sine)',...
	%			4 = 	 'fixed response   (Box-car)'};

% epoch-related stuff

TOR.epoch.resptype = 4;
	% 1 cosine
	% 2 exponential
	% 3 half-sine
	% 4 box-car

TOR.epoch.convhrf = 1;
	% convolve with hrf: yes or no, 1 or 0.

TOR.epoch.addtd = 0;
	% add temporal derivatives: yes or no, 1 or 0.

TOR.epoch.elength = [16 16];
	% epoch length in scans for each trial type
	
TOR.c=c ;
TOR.v = length(c) ./ TOR.numscans;		
	% number of conditions/trial types in each scan

TOR.Cname = {'SCSLL' 'NCSLL' 'SCNLL' 'NCNLL' 'SCSLR' 'NCSLR' 'SCNLR' 'NCNLR' 'probe1' 'probe2' 'firstL' 'firstR'};	
% Switch Counter No-switch Location Right vis field = SCNLR
	% names for each trial type

TOR.Volterra = 0;
TOR.userc = 0;
TOR.Global = 'Scaling';
TOR.cLF = 'specify';									
												% high-pass filter:
												% 'specify' or 'none'
												
TOR.HParam = [80 80 80 80 80 80 80 80];
	% high pass filter cutoff for each session

TOR.LFstr = num2str(TOR.HParam);						% high-pass filter text string - in bottom of design window
TOR.cHF = 'hrf';										% low-pass filter option
TOR.cVi = 'none';										% none or 'AR(1)'
TOR.bFcon = 0;
TOR.bEstNow = 1;	
	% do not estimate, so that we may modify CFG.mat to insert bigmask!
	% but this is changed automatically below, so no need to worry if TOR.dobigmask = 1;
	% estimate now?

% User-specified regressors (realignment parameters)
% TOR.usregs = real';
% TOR.usna = {'x  ','y  ','z  ','rol','pit','yaw'};

% Kalina Christoff's bigmask stuff for modifying analysis mask
% -------------------------------------------------------------------
TOR.dobigmask = 1;			% insert bigmask before estimating model
TOR.in_create = 'y';		% create bigmask image file
TOR.in_segment = 'n';		% use already created segmented images if 'n'
TOR.in_modify = 'y';		% modify spmCFG.mat before estimating
TOR.in_estim = 'y';			% estimate after modifying spmCFG.mat
TOR.PG = EXPT.SUBJECT.canonicalT1;
	% template image (T1.img or T2.img)
TOR.PF = EXPT.SUBJECT.nspgr(find(strcmp(EXPT.SUBJECT.subdir,SubjCode)),:);				
	% image to segment
TOR.cfg_file = [TOR.resdir '/SPMcfg.mat'];
		% cfg.mat file to modify
TOR.bigmaskimg = [EXPT.SUBJECT.studydir filesep subs{1} '/anatomy/bigmask.img'];

   clear a, clear b, clear c,clear contrast, clear names
  
% define contrast matrices
% ====================================================================== 

%	{'SCSLL' 'NCSLL' 'SCNLL' 'NCNLL' 'SCSLR' 'NCSLR' 'SCNLR' 'NCNLR' 'probe1' 'probe2' 'firstL' 'firstR'};	

cnum = 1;
a = [1 1 1 1 1 1 1 1 0 0 0 0];
c{cnum} = [a a a a a a a a];
names{cnum} = 'All effects';
type{cnum} = 'F';

cnum = 2; % 0003
a = [1 1 1 1 -1 -1 -1 -1 0 0 0 0];
c{cnum} = [a a a a a a a a];
names{cnum} = 'L vs R';
type{cnum} = 'T';

cnum = 3;
a = [-1 -1 -1 -1 1 1 1 1 0 0 0 0];
c{cnum} = [a a a a a a a a];
names{cnum} = 'R vs L';
type{cnum} = 'T';

cnum = 4; % 0005
a = [1 -1 1 -1 1 -1 1 -1 0 0 0 0];
c{cnum} = [a a a a a a a a];
names{cnum} = 'Counter Switch Act';
type{cnum} = 'T';

cnum = 5;
a = [-1 1 -1 1 -1 1 -1 1 0 0 0 0];
c{cnum} = [a a a a a a a a];
names{cnum} = 'Counter Switch Deact';
type{cnum} = 'T';

cnum = 6;
a = [-1 1 -1 1 -1 1 -1 1 0 0 0 0];
c{cnum} = [a a a a a a a a];
names{cnum} = 'Counter Switch Deact';
type{cnum} = 'T';

cnum = 7; % 0008
a = [1 1 -1 -1 1 1 -1 -1 0 0 0 0];
c{cnum} = [a a a a a a a a];
names{cnum} = 'Loc Switch Act';
type{cnum} = 'T';

cnum = 8;
a = [-1 -1 1 1 -1 -1 1 1 0 0 0 0];
c{cnum} = [a a a a a a a a];
names{cnum} = 'Loc Switch Deact';
type{cnum} = 'T';

cnum = 9; 
a = [1 1 1 1 -1 -1 -1 -1 0 0 0 0];
b = [-1 -1 -1 -1 1 1 1 1 0 0 0 0];
c{cnum} = [a a a a b b b b];
names{cnum} = 'L vs R by half';
type{cnum} = 'T';

cnum = 10; 
a = [1 1 1 1 -1 -1 -1 -1 0 0 0 0];
b = [-1 -1 -1 -1 1 1 1 1 0 0 0 0];
c{cnum} = [b b b b a a a a ];
names{cnum} = 'R vs L by half';
type{cnum} = 'T';

cnum = 11; 
a = [1 -1 1 -1 1 -1 1 -1 0 0 0 0];
b = -a
c{cnum} = [a a a a b b b b];
names{cnum} = 'Counter Sw by Half';
type{cnum} = 'T';

cnum = 12; 
a = [1 1 -1 -1 1 1 -1 -1 0 0 0 0];
b = -a
c{cnum} = [a a a a b b b b];
names{cnum} = 'Loc Sw by Half';
type{cnum} = 'T';


% define results thresholding choices
% [run tor_spm_results_ui(TOR)]
% ====================================================================== 

TOR.maskWOtherContrasts = 0;
	% mask with other contrasts; 1 or 0

TOR.multCompCorrect = 'uncorrected'; 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

TOR.u = .005;
	% height threshold - T or p value

TOR.k = 10;
	% extent threshold - number of voxels

TOR.Ic = 'all';
	% optional: which contrast to show in results (number)
	% or 'all'.  If 'all', uses all contrasts, in order.
	% if using SelectedContrast, below, this is automatically changed.

% also uses TOR.resdir


% For estimating only a single contrast
% ======================================================================
SelectedContrast = 0;	% estimate only a single contrast if 1, otherwise all
whichContrast = 34;		% number of contrast to estimate and get results for




% ====================================================================== 
% * END USER INPUT
% ====================================================================== 


EXPT.SUBJECT.contrasts = c;
EXPT.SUBJECT.connames = names;
EXPT.SUBJECT.ctype = type;
EXPT.SUBJECT.STATINPUT = TOR;



% For estimating only a single contrast
% ======================================================================
if SelectedContrast
	myContrastName = names{whichContrast};
	TOR.Ic = whichContrast;
else
    myContrastName = [];
end

% make sure we're in the correct directory
% ====================================================================== 
TOR
if ~(exist(TOR.resdir) == 7),
	disp(['Making directory ' TOR.resdir])
	eval(['!mkdir ' TOR.resdir])
end
eval(['cd ' TOR.resdir])
cont = [];
if ~(exist('goOK') == 1), goOK = 0; end
if ~(goOK) == 1
	while ~(strcmp(cont,'y') | strcmp('cont','n'))
		cont = input(['Data will be saved in ' pwd '. Save results here (y\n)?'],'s');
		if strcmp(cont,'n'),error('Please change to the correct results directory.'),end
	end
end


% get the list of actual file names
% ====================================================================== 

for i = 1:length(TOR.imgdir)
	[fNames,dummy] = spm_list_files(TOR.imgdir{i},TOR.imgwildcard);
	% ...and add the directory name
 	a = repmat([TOR.imgdir{i} filesep],length(fNames),1);
	fNames = [a fNames];
	TOR.imgnames{i} = fNames;
	if isempty(fNames),error(['No files matching ' TOR.imgwildcard ' in ' TOR.imgdir{i}]),end 
end
 
% save parameter structure as .mat file
% ====================================================================== 
save Stats_Input_Params TOR   

   

% set up bigmask insertion stuff
% ====================================================================== 
if TOR.dobigmask
	TOR.bEstNow = 0;
end



% build and estimate the model (or just set up)
% ====================================================================== 
if TOR.estimatemodel
	tor_spm_fmri_spm_ui(TOR)
end  


% do CFG modification for bigmask, and estimate if specified
% ====================================================================== 
if TOR.dobigmask & TOR.estimatemodel
	tor_glm_specmask(TOR)
end


% save parameter structure and contrast structure as .mat file
% ====================================================================== 
save Stats_Input_Params TOR c names type  
  



% get the contrasts and save in xCon.mat  
% ==============================================================================
if TOR.estimatecontrasts
	estimateContrasts(c,names,type,SelectedContrast)
end


% threshold and display results
% ====================================================================== 

% top part is for getting results of ALL contrasts
% bottom part (last line) is for getting results and choosing a contrast w/ the GUI
%	or choosing the single contrast number specified in TOR.Ic

if TOR.getresults

if isfield(TOR,'Ic')
	if strcmp(TOR.Ic,'all')
		load xCon
		for i = 1:length(xCon)
			TOR.Ic = i;
			[hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(TOR);
			spm_list('List',SPM,VOL,[],[],'',hReg);				% list volume stats in window
			spm_print
			pause(5)
		end

		eval(['!mv spm99.ps ' SubjCode '_results.ps'])
	end
else

	[hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(TOR);
	spm_list('List',SPM,VOL,[],[],'',hReg);							% list volume stats in window
end

end


diary off

