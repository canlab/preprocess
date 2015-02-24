function [xX,Sess] = tor_spm_fmri_spm_ui(TOR)
% Setting up the general linear model for fMRI time-series
% FORMAT [xX,Sess] = spm_fmri_spm_ui
%
% xX            - structure describing design matrix
% xX.X          - design matrix
% xX.dt         - time bin {secs}
% xX.TOR.RT         - Repetition time {secs}
% xX.iH         - vector of H partition (condition effects)      indices,
% xX.iC         - vector of C partition (covariates of interest) indices
% xX.iB         - vector of B partition (block effects)          indices
% xX.iG         - vector of G partition (nuisance variables)     indices
% xX.Xnames     - cellstr of effect names corresponding to columns
%                 of the design matrix
%
% Sess{s}.BFstr   - basis function description string
% Sess{s}.DSstr   - Design description string
% Sess{s}.row     - scan   indices      for session s
% Sess{s}.col     - effect indices      for session s
% Sess{s}.name{i} - of ith trial type   for session s
% Sess{s}.ind{i}  - column indices      for ith trial type {within session}
% Sess{s}.bf{i}   - basis functions     for ith trial type
% Sess{s}.sf{i}   - stick functions     for ith trial type
% Sess{s}.ons{i}  - stimuli onset times for ith trial type (secs)
% Sess{s}.pst{i}  - peristimulus times  for ith trial type (secs)
% Sess{s}.para{i} - vector of paramters for ith trial type
%____________________________________________________________________________
%
% spm_fmri_spm_ui configures the design matrix, data specification and
% filtering that specify the ensuing statistical analysis. These
% arguments are passed to spm_spm that then performs the actual parameter
% estimation.
%
% The design matrix defines the experimental design and the nature of
% hypothesis testing to be implemented.  The design matrix has one row
% for each scan and one column for each effect or explanatory variable.
% (e.g. regressor or stimulus function).  The parameters are estimated in
% a least squares sense using the general linear model.  Specific profiles
% within these parameters are tested using a linear compound or contrast
% with the T or F statistic.  The resulting statistical map constitutes 
% an SPM.  The SPM{T}/{F} is then characterized in terms of focal or regional
% differences by assuming that (under the null hypothesis) the components of
% the SPM (i.e. residual fields) behave as smooth stationary Gaussian fields.
%
% spm_fmri_spm_ui allows you to (i) specify a statistical model in terms
% of a design matrix, (ii) review that design, (iii) associate some data
% with a pre-specified design [or (iv) specify both the data and design]
% and then proceed to estimate the parameters of the model.
% Inferences can be made about the ensuing parameter estimates (at a first
% or fixed-effect level) in the results section, or they can be re-entered
% into a second (random-effect) level analysis by treating the session or 
% subject-specific [contrasts of] parameter estimates as new summary data.
% Inferences at any level obtain by specifying appropriate T or F contrasts
% in the results section to produce SPMs and tables of p values and statistics.
%
% spm_fmri_spm calls spm_fMRI_design which allows you to configure a
% design matrix in terms of events or epochs.  This design matrix can be
% specified before or during data specification.  In some instances
% (e.g. with stochastic designs that have to realized before data
% acquisition) it is necessary to build the design matrix first and then
% select the corresponding data.  In others it may be simpler to specify
% the data and then the design.  Both options are supported.  Once the
% design matrix, data and filtering have been specified spm_fmri_spm_ui
% calls spm_spm to estimate the model parameters that are then saved for
% subsequent analysis.
%
% spm_fMRI_design allows you to build design matrices with separable
% session-specific partitions.  Each partition may be the same (in which
% case it is only necessary to specify it once) or different.  Responses
% can be either event- or epoch related, where the latter model prolonged
% and possibly time-varying responses to state-related changes in
% experimental conditions.  Event-related response are modelled in terms
% of responses to instantaneous events.  Mathematically they are both
% modeled by convolving a series of delta (or stick) functions,
% indicating the onset of an event or epoch with a set of basis
% functions.  These basis functions can be very simple, like a box car,
% or may model voxel-specific forms of evoked responses with a linear
% combination of several basis functions (e.g. a Fourier set).  Basis
% functions can be used to plot estimated responses to single events or
% epochs once the parameters (i.e. basis function coefficients) have
% been estimated.  The importance of basis functions is that they provide
% a graceful transition between simple fixed response models (like the
% box-car) and finite impulse response (FIR) models, where there is one
% basis function for each scan following an event or epoch onset.  The
% nice thing about basis functions, compared to FIR models, is that data
% sampling and stimulus presentation does not have to be sychronized
% thereby allowing a uniform and unbiased sampling of peri-stimulus time.
% 
% Event-related designs may be stochastic or deterministic.  Stochastic
% designs involve one of a number of trial-types occurring with a
% specified probably at successive intervals in time.  These
% probabilities can be fixed (stationary designs) or time-dependent
% (modulated or non-stationary designs).  The most efficient designs
% obtain when the probabilities of every trial type are equal and this is
% enforced in SPM.  The modulation of non-stationary designs is simply
% sinusoidal with a period of 32 seconds.  A critical aspect of
% stochastic event-related designs is whether to include null events or
% not.  If you wish to estimate the evoke response to a specific event
% type (as opposed to differential responses) then a null event must be
% included (even though it is not modelled explicitly).
% 
% The choice of basis functions depends upon the nature of the inference
% sought.  One important consideration is whether you want to make
% inferences about compounds of parameters (i.e.  contrasts).  This is
% the case if (i) you wish to use a SPM{T} to look separately at
% activations and deactivations or (ii) you with to proceed to a second
% (random-effect) level of analysis.  If this is the case then (for
% event-related studies) use a canonical hemodynamic response function
% (HRF) and derivatives with respect to latency (and dispersion).  Unlike
% other bases, contrasts of these effects have a physical interpretation
% and represent a parsimonious way of characterising event-related
% responses.  Bases such as a Fourier set require the SPM{F} for
% inference and preclude second level analyses.
% 
% See spm_fMRI_design for more details about how designs are specified.
%
% Serial correlations in fast fMRI time-series are dealt with as
% described in spm_spm.  At this stage you need to specific the filtering
% that will be applied to the data (and design matrix).  This filtering
% is important to ensure that bias in estimates of the standard error are
% minimized.  This bias results from a discrepancy between the estimated
% (or assumed) auto-correlation structure of the data and the actual
% intrinsic correlations.  The intrinsic correlations will be estimated
% automatically using an AR(1) model during parameter estimation.  The
% discrepancy between estimated and actual intrinsic (i.e. prior to
% filtering) correlations are greatest at low frequencies.  Therefore
% specification of the high-pass component of the filter is particularly
% important.  High pass filtering is now implemented at the level of the
% filtering matrix K (as opposed to entering as confounds in the design
% matrix).  The default cutoff period is twice the maximum time interval
% between the most frequently occurring event or epoch (i.e the minium of
% all maximum intervals over event or epochs).
%
% Burst Mode is a specialist design for intermittent epochs of acquisitions
% (used for example to allow for intercalated EEG recording).  Each burst
% is treated as a session but consistent within-session effects (e.g. T1
% effects) are modeled in X.bX.  The primary use of this mode is to generate
% parameter estimate images for a second level analysis.
%
% N.B: Burst Mode is an unsupported 'in-house' developmental feature.
%      The feature is currently hidden from the user.
%
%-----------------------------------------------------------------------
% Refs:
%
% Friston KJ, Holmes A, Poline J-B, Grasby PJ, Williams SCR, Frackowiak
% RSJ & Turner R (1995) Analysis of fMRI time-series revisited. NeuroImage
% 2:45-53
%
% Worsley KJ and Friston KJ (1995) Analysis of fMRI time-series revisited -
% again. NeuroImage 2:178-181
%
% Friston KJ, Frith CD, Frackowiak RSJ, & Turner R (1995) Characterising
% dynamic brain responses with fMRI: A multivariate approach NeuroImage -
% 2:166-172
%
% Frith CD, Turner R & Frackowiak RSJ (1995) Characterising evoked 
% hemodynamics with fMRI Friston KJ, NeuroImage 2:157-165
%
% Josephs O, Turner R and Friston KJ (1997) Event-related fMRI, Hum. Brain
% Map. 0:00-00
%
%_______________________________________________________________________
% @(#)spm_fmri_spm_ui.m	2.31 Karl Friston, Jean-Baptiste Poline, Christian Buchel 99/10/14

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'model',index_of_Analysis};
%_______________________________________________________________________

global torindex

SCCSid  = '2.31';

global BCH; %- used as a flag to know if we are in batch mode or not.

%-GUI setup
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','fMRI stats model setup',0);
spm_help('!ContextHelp',mfilename)


% get design matrix and/or data
%=======================================================================
MType = {'specify a model',...
	   'review a specified model',...
	   'estimate a specified model',...
	   'specify and estimate a model'};
%TOR.MT    = spm_input('What would you like to do?',1,'m',MType,...
%                  'batch',{},'types');
% tor               

%-Initialise output arguments in case return early
xX   = [];
Sess = [];

switch TOR.MT
%-----------------------------------------------------------------------

	case 1
	% specify a design matrix
	%---------------------------------------------------------------
	if sf_abort, spm_clf(Finter), return, end
	[xX,Sess] = spm_fMRI_design;
	spm_fMRI_design_show(xX,Sess);
	return

	case 2
	% 'review a specified model'
	%---------------------------------------------------------------
	spm_clf(Finter)
	[xX,Sess] = spm_fMRI_design_show;
	return

	case 3
	% load pre-specified design matrix
	%---------------------------------------------------------------
	if sf_abort([2,3]), spm_clf(Finter), return, end
	if isempty(BCH)
	   % load(spm_get(1,'fMRIDesMtx.mat','Select SPM_fMRIDesMtx.mat'));
	   % Tor: automatically load one in current directory
	   load('SPM_fMRIDesMtx.mat');
	else
	   load('SPM_fMRIDesMtx.mat');
	end


	% get filenames
	%---------------------------------------------------------------
	nsess  = length(xX.iB);
	nscan  = zeros(1,nsess);
	for  i = 1:nsess
		nscan(i) = length(find(xX.X(:,xX.iB(i))));
	end
	P      = [];
	if nsess < 16
		for i = 1:nsess
		   str = sprintf('select scans for session %0.0f',i);
		   if isempty(BCH)
			% q = spm_get(Inf,'.img',str);
			% Tor: read from input structure
			q = TOR.imgnames{i};
		   else
			q = sf_bch_get_q(i);
		   end %- 
         P   = strvcat(P,q);
         
		end
	else
		str   = sprintf('select scans for this study');
		if isempty(BCH)
			P     = spm_get(sum(nscan),'.img',str);
		else
		   for i = 1:nsess
			q = sf_bch_get_q(i);
			P = strvcat(P,q);
		   end
		end
	end
   
   % tor               
   P(1,:)
   
   
	% Repeat time
	%---------------------------------------------------------------
	% This should be ok without this - uses TOR.RT as input to design matrix maker.
	%TOR.RT     = xX.TOR.RT;
	%TOR.RT	

	case 4
	% get filenames and design matrix
	%---------------------------------------------------------------
	if sf_abort, spm_clf(Finter), return, end
	spm_input('Scans & sessions...',1,'d',mfilename,'batch')
	%nsess  = spm_input(['number of sessions'],'+1','e',1,...
   %   'batch',{},'nsess');
   
   % tor               
   nsess = TOR.numscans;
   
	nscan  = zeros(1,nsess);
	P      = [];
	for  i = 1:nsess
      		str  = sprintf('select scans for session %0.0f',i);
            
      		if isempty(BCH)
         		% tor
      			% r = getfiles([TOR.imgdir filesep 'scan' num2str(i) filesep TOR.imgnames]);
			% old - looks in scan subdirs automatically.

			q = TOR.imgnames{i};
	
			% new: input is actual list of image names.
            		% don't need this to convert - old stuff.
            		%r{1}
            		%for z = 1:size(r,1)
            		%   q(z,:) = r{z};
            		%end
            
		   %q = spm_get(Inf,'.img',str);
		else
		   q = sf_bch_get_q(i);
		end
 		P        = strvcat(P,q);
		nscan(i) = size(q,1);
	end
   
   % tor               
   P(1,:)
   
	% get Repeat time
	%---------------------------------------------------------------
	% TOR.RT  = spm_input('Interscan interval {secs}','+1','batch',{},'RT');
   	% defined earlier.
   % tor               
   
   
	% get design matrix
	%---------------------------------------------------------------
	[xX,Sess] = tor_spm_fMRI_design(nscan,TOR.RT,TOR);

end

% Assemble other deisgn parameters
%=======================================================================
spm_help('!ContextHelp',mfilename)
spm_input('Global intensity normalisation...',1,'d',mfilename,'batch')

% get rows
%-----------------------------------------------------------------------
for i = 1:nsess
	row{i} = find(xX.X(:,xX.iB(i)));
end
BFstr  = Sess{1}.BFstr;
DSstr  = Sess{1}.DSstr;


% Global normalization
%-----------------------------------------------------------------------
str    = 'remove Global effects';
%Global = spm_input(str,'+1','scale|none',{'Scaling' 'None'},...
%   'batch',{},'global_effects');
Global = TOR.Global;
if ischar(Global),
	Global = {Global};
end


% Burst mode
%-----------------------------------------------------------------------
% if length(nscan) > 16 & ~any(diff(nscan))
% 	spm_input('Burst mode?','+1','d',mfilename,'batch')
% 	BM    = spm_input('Burst mode?','+1','y/n',[1 0],2,...
% 			    'batch',{},'burst_mode');
% else
% 	BM    = 0;
% end

%-** Burst mode is a developmental feature, disabled due to problems!
%-** See (KJF) http://www.mailbase.ac.uk/lists/spm/1999-10/0009.html
BM = 0;

% Model scan effects as oppsed to session effects if burst mode
%-----------------------------------------------------------------------
if BM
	k         = nscan(1);
	l         = length(xX.iC);
	xX.X      = xX.X(:,xX.iC);
	xX.Xnames = xX.Xnames(xX.iC);
	xX.X      = [xX.X kron(ones(nsess,1),eye(k))];
	xX.iB     = [1:k] + l;
	for   i = 1:k
		X.Xnames{l + i} = sprintf('scan: %i ',i);
	end
	DSstr = '[burst-mode]';
end


% Temporal filtering
%=======================================================================
spm_input('Temporal autocorrelation options','+1','d',mfilename,'batch')

% High-pass filtering
%-----------------------------------------------------------------------
if BM
	cLF = 'none';
else
	%cLF = spm_input('High-pass filter?','+1','b','none|specify',...
   %		'batch',{},'HF_fil');
   cLF = TOR.cLF;
end
   
switch cLF

	case 'specify'

	% default based on peristimulus time
	% param = cut-off period (max = 512, min = 32)
	%---------------------------------------------------------------
	HParam = 512*ones(1,nsess);
	for  i = 1:nsess
		for j = 1:length(Sess{i}.pst)
		   HParam(i) = min([HParam(i) 2*max(TOR.RT + Sess{i}.pst{j})]);
		end
	end
	HParam = ceil(HParam);
	HParam(HParam < 32) = 32;
	str    = 'session cutoff period (secs)';
	%HParam = spm_input(str,'+1','e',HParam,[1 nsess],...
	%                  'batch',{},'HF_cut');
   HParam = TOR.HParam;
   
	% LF description
	%---------------------------------------------------------------
	%LFstr = sprintf('[min] Cutoff period %d seconds',min(HParam));
   LFstr = TOR.LFstr;
   
	case 'none'
	%---------------------------------------------------------------
	HParam = cell(1,nsess);
	LFstr  = cLF;

end


% Low-pass filtering
%-----------------------------------------------------------------------
if BM
	cHF = 'none';
else
	%cHF = spm_input('Low-pass filter?','+1','none|Gaussian|hrf',...
	%		'batch',{},'LF_fil');
cHF = TOR.cHF;

end
switch cHF

	case 'Gaussian'
	%---------------------------------------------------------------
	LParam  = spm_input('Gaussian FWHM (secs)','+1','r',4,...
			    'batch',{},'LF_cut');
	HFstr   = sprintf('Gaussian FWHM %0.1f seconds',LParam);
	LParam  = LParam/sqrt(8*log(2));

	case {'hrf', 'none'}
	%---------------------------------------------------------------
	LParam  = [];
	HFstr   = cHF;

end

% create filter struct and band-pass specification
%-----------------------------------------------------------------------
for i = 1:nsess
	K{i} = struct(	'HChoice',	cLF,...
			'HParam',	HParam(i),...
			'LChoice',	cHF,...
			'LParam',	LParam,...
			'row',		row{i},...
			'RT',		TOR.RT);
end


% intrinsic autocorrelations (Vi)
%-----------------------------------------------------------------------
str     = 'Model intrinsic correlations?';
cVimenu = {'none','AR(1)'};
%cVi     = spm_input(str,'+1','b',cVimenu,'batch',{},'int_corr');
cVi = TOR.cVi;

%-Estimation options
%=======================================================================
spm_input('Estimation options',1,'d',mfilename,'batch')

%-Generate default trial-specific F-contrasts specified by session?
%-----------------------------------------------------------------------
%bFcon = spm_input('Setup trial-specific F-contrasts?','+1','y/n',[1,0],1,...
%		'batch',{},'trial_fcon');
bFcon = TOR.bFcon;

%-Estimate now or later?
%-----------------------------------------------------------------------
%bEstNow = spm_input('estimate?','_','b','now|later',[1,0],1,...
%		'batch',{},'now_later');
bEstNow = TOR.bEstNow;

%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spm_clf(Finter);
spm('FigName','Configuring, please wait...',Finter,CmdLine);
spm('Pointer','Watch');


% Contruct K and Vi structs
%=======================================================================
K       = spm_filter('set',K);


% Adjust for missing scans
%-----------------------------------------------------------------------
[xX,Sess,K,P,nscan,row] = spm_bch_tsampl(xX,Sess,K,P,nscan,row); %-SR

% create Vi struct
%-----------------------------------------------------------------------
if isfield(TOR,'Vi'),
   Vi = TOR.Vi;
else
   %spm's way: identity matrix
   Vi      = speye(sum(nscan));
end
xVi     = struct('Vi',Vi,'Form',cVi);
% according to the help, including row makes it use AR(1) model
% so I'm leaving it out.
%for   i = 1:nsess
%	xVi.row{i} = row{i};
%end

% get file identifiers and Global values
%=======================================================================
fprintf('%-40s: ','Mapping files')                                   %-#
VY     = spm_vol(P);
fprintf('%30s\n','...done')                                          %-#

if any(any(diff(cat(1,VY.dim),1,1),1)&[1,1,1,0])
	error('images do not all have the same dimensions'),           end
if any(any(any(diff(cat(3,VY.mat),1,3),3)))
	error('images do not all have same orientation & voxel size'), end


%-Compute Global variate
%-----------------------------------------------------------------------
GM     = 100;
q      = sum(nscan);
g      = zeros(q,1);
fprintf('%-40s: %30s','Calculating globals',' ')                     %-#
for i  = 1:q
    disp(int2str(i))
	fprintf('%s%30s',sprintf('\b')*ones(1,30),sprintf('%4d/%-4d',i,q)) %-#
	g(i) = spm_global(VY(i));
end
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


% scale if specified (otherwise session specific grand mean scaling)
%-----------------------------------------------------------------------
gSF     = GM./g;
if strcmp(Global,'None')
	for i = 1:nsess
		j      = row{i};
		gSF(j) = GM./mean(g(j));
	end
end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%-----------------------------------------------------------------------
for  i = 1:q, VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*gSF(i); end


%-Masking structure
%-----------------------------------------------------------------------
xM     = struct('T',	ones(q,1),...
		'TH',	g.*gSF,...
		'I',	0,...
		'VM',	{[]},...
		'xs',	struct('Masking','analysis threshold'));


%-Complete design matrix (xX)
%=======================================================================
xX.K   = K;
xX.xVi = xVi;



%-Effects designated "of interest" - constuct F-contrast structure array
%-----------------------------------------------------------------------
if length(xX.iC)
	F_iX0  = struct(	'iX0',		xX.iB,...
				'name',		'effects of interest');
else
	F_iX0  = [];
	DSstr  = 'Block [session] effects only';
end

%-Trial-specifc effects specified by Sess
%-----------------------------------------------------------------------
%-NB: With many sessions, these default F-contrasts can make xCon huge!
if bFcon
	i      = length(F_iX0) + 1;
	if (Sess{1}.rep)
		for t = 1:length(Sess{1}.name)
			u     = [];
			for s = 1:length(Sess)
				u = [u Sess{s}.col(Sess{s}.ind{t})];
			end
			q             = 1:size(xX.X,2);
			q(u)          = [];
			F_iX0(i).iX0  = q;
			F_iX0(i).name = Sess{s}.name{t};
			i             = i + 1;
		end
	else
		for s = 1:length(Sess)
			str   = sprintf('Session %d: ',s);
			for t = 1:length(Sess{s}.name)
				q             = 1:size(xX.X,2);
				q(Sess{s}.col(Sess{s}.ind{t})) = [];
				F_iX0(i).iX0  = q;
				F_iX0(i).name = [str Sess{s}.name{t}];
				i             = i + 1;
			end
		end
	end
end


%-Design description (an nx2 cellstr) - for saving and display
%=======================================================================
for i    = 1:length(Sess), ntr(i) = length(Sess{i}.name); end
sGXcalc  = 'mean voxel value';
sGMsca   = 'session specific';
xsDes    = struct(	'Design',			DSstr,...
			'Basis_functions',		BFstr,...
			'Number_of_sessions',		sprintf('%d',nsess),...
			'Conditions_per_session',	sprintf('%-3d',ntr),...
			'Interscan_interval',		sprintf('%0.2f',TOR.RT),...
			'High_pass_Filter',		LFstr,...
			'Low_pass_Filter',		HFstr,...
			'Intrinsic_correlations',	xVi.Form,...
			'Global_calculation',		sGXcalc,...
			'Grand_mean_scaling',		sGMsca,...
			'Global_normalisation',		Global);
%-global structure
%-----------------------------------------------------------------------
xGX.iGXcalc  = Global{1};
xGX.sGXcalc  = sGXcalc;
xGX.rg       = g;
xGX.sGMsca   = sGMsca;
xGX.GM       = GM;
xGX.gSF      = gSF;


%-Save SPMcfg.mat file
%-----------------------------------------------------------------------
fprintf('%-40s: ','Saving SPMstats configuration')                   %-#
save SPMcfg SPMid xsDes VY xX xM xGX F_iX0 Sess
fprintf('%30s\n','...SPMcfg.mat saved')                              %-#

 
%-Display Design report
%=======================================================================
fprintf('%-40s: ','Design reporting')                                %-#
spm_DesRep('DesMtx',xX,reshape({VY.fname},size(VY)),xsDes)
fprintf('%30s\n','...done')                                          %-#



%-Analysis Proper
%=======================================================================
spm_clf(Finter);
spm('FigName','fMRI stats models',Finter,CmdLine);
if bEstNow
	spm('Pointer','Watch')
	spm('FigName','Stats: estimating...',Finter,CmdLine);
	spm_spm(VY,xX,xM,F_iX0,Sess,xsDes);
	spm('Pointer','Arrow')
else
	spm_clf(Finter)
	spm('FigName','Stats: configured',Finter,CmdLine);
	spm('Pointer','Arrow')
	spm_DesRep('DesRepUI',struct(	'xX',		xX,...
					'VY',		VY,...
					'xM',		xM,...
					'F_iX0',	F_iX0,...
					'Sess',		{Sess},...
					'xsDes',	xsDes,...
					'swd',		pwd,...
					'SPMid',	SPMid,...
					'cfg',		'SPMcfg'));
end


%-End: Cleanup GUI
%-----------------------------------------------------------------------
fprintf('\n\n')



%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function abort = sf_abort(i)
%=======================================================================

if nargin<1, i=[1:3]; end
tmp    = zeros(1,3);
tmp(i) = 1;
tmp = tmp & [	exist(fullfile('.','SPM_fMRIDesMtx.mat'),'file')==2 ,...
		exist(fullfile('.','SPMcfg.mat'),        'file')==2 ,...
		exist(fullfile('.','SPM.mat'),           'file')==2 ];
if any(tmp)
	str = {	'    SPM fMRI design matrix definition (SPM_fMRIDesMtx.mat)',...
		'    SPMstats configuration            (SPMcfg.mat)',...
		'    SPMstats results files            (inc. SPM.mat)'};
	str = {	'Current directory contains existing SPMstats files:',...
		str{tmp},['(pwd = ',pwd,')'],' ',...
		'Continuing will overwrite existing files!'};

	abort = spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename,...
                         'batch',{},'stop_writing');
	if abort, fprintf('%-40s: %30s\n\n',...
		'Abort...   (existing SPMstats files)',spm('time')), end
else
	abort = 0;
end

function q = sf_bch_get_q(i)
%=======================================================================
% This is to deal with a specific case where the sampling   
% isn't regular. Only implemented in bch mode.
% 
q 	= spm_input('batch',{},'files',i);
files 	= q;
t_sampl  = spm_input('batch',{},'time_sampl',i);
remain 	= spm_input('batch',{},'remain',i);
q(remain,:) = files;

%- fills the gap with the first images .... 
%- there should be enough first images !
q(t_sampl,:) = files(1:length(t_sampl),:);
