% global fmriTEST fmriDIR;             %set global variabels 
% fmriTEST = 0; 
% fmriDIR  = '/data4/intext2';
% f99_CWD('LOG');                   %make the LOG and RESULTS directory 


% make sure we're in the correct directory
% ====================================================================== 
cont = [];
while ~(strcmp(cont,'y') | strcmp('cont','n'))
	cont = input(['Current directory is ' pwd '. Save results here (y\n)?'],'s');
	if strcmp(cont,'n'),error('Please change to the correct results directory.'),end
end

myDayCode = input('Enter day code ','s');

% open a diary file to save output
% ====================================================================== 
diary Stat_Script_Diary

% copy the design matrix down to this directory
% ====================================================================== 
!cp /data/ektrain/DesMtxGaravan.mat ./SPM_fMRIDesMtx.mat


% define the vectors of onsets for each trial type
% ====================================================================== 

% For Garavan
%	Baseline	HiSwitch	LoSwitch
c{1} = [0 340];
c{2} = [20:40:300];
c{3} = [40:40:320];



% define input parameters to spm_spm_ui
% ====================================================================== 

TOR.refslice = 1;
TOR.fmriT = 16;
TOR.RT=1;

TOR	.MT=3;												
	% 1 = specify
	% 2 = review
	% 3 = estimate
	% 4 = specify and estimate


TOR.imgdir = ['/data/ektrain/task/' myDayCode '/garavanarrow'];
TOR.resdir = [pwd];
TOR.imgwildcard = 'srra*img';
  

TOR.vardurations=0;   
TOR.Ptype='none'; 
TOR.Etype='linear';   
TOR.Rov='epochs';

TOR.Cov=4;
	% For Events:	1 = hrf alone, 2 = hrf + td
	% For Epochs:	1 = 'basis functions  (Discrete Cosine Set)',...
	%			2 = 'basis functions  (Mean & exponential decay)',...
	%			3 = 	 'fixed response   (Half-sine)',...
	%			4 = 	 'fixed response   (Box-car)'};
	
TOR.c=c ;
TOR.v = 3;		
	% number of conditions/trial types in each scan

TOR.Cname = {'Baseline','HiSwitch','LoSwitch'};	
	% names for each trial type

TOR.Volterra = 0;
TOR.userc = 0;
TOR.Global = 'None';
TOR.cLF = 'specify';
TOR.HParam = [82];
	% high pass filter cutoff for each session

TOR.LFstr = 'none';										% high-pass filter text string - not used? - just a text message...
TOR.cHF = 'hrf';										% low-pass filter option
TOR.cVi = 'none';										% none or 'AR(1)'
TOR.bFcon = 0;
TOR.bEstNow = 1;
	% estimate now?

% User-specified regressors (realignment parameters)
% TOR.usregs = real';
% TOR.usna = {'x  ','y  ','z  ','rol','pit','yaw'};


% get the list of actual file names
% ====================================================================== 

[TOR.imgnames,dummy] = spm_list_files(TOR.imgdir,TOR.imgwildcard);
% ...and add the directory name
 a = repmat([TOR.imgdir filesep],length(TOR.imgnames),1);
TOR.imgnames = [a TOR.imgnames];

 
% save parameter structure as .mat file
% ====================================================================== 
save Stats_Input_Params TOR   

   
% build and estimate the model
% ====================================================================== 
tor_spm_fmri_spm_ui(TOR)
   



   clear a, clear b, clear c,clear contrast, clear names
   

% define contrast matrices
% ==============================================================================
  c{1} = [0 1 -1];
  names{1} = 'HiSwitch - LoSwitch';
  type{1} = 'T';

  c{2} = [0 -1 1];
  names{2} = 'LoSwitch - HiSwitch';
  type{2} = 'T';

  c{3} = [-1 1 0];
  names{3} = 'HiSwitch - Baseline';
  type{3} = 'T';

  c{4} = [-1 0 1];
  names{4} = 'LoSwitch - Baseline';
  type{4} = 'T';

  c{5} = [-2 1 1];
  names{5} = 'BOSwitch - Baseline';
  type{5} = 'T';

  c{6} = [2 -1 -1];
  names{6} = 'Baseline - BOSwitch';
  type{6} = 'T';


% save parameter structure and contrast structure as .mat file
% ====================================================================== 
save Stats_Input_Params TOR c names type  
  

 % get the contrasts and save in xCon.mat  
 % ==============================================================================

 
estimateContrasts

diary off
