function f99_spm_fmri_calc(DIRe,EXP,kulRES,glob,BM,cLF,HParam_in,cHF,LParam_in,cVi,bFcon)

% DIRe           - [matrix of directories]
% EXP            - [matrix of experiment names]
% RES            - directory where to find SPM_fMRIDesMtx
%                  directory where to write the results
% glob           - global normalisation ('Scaling' or 'None')
% BM             - burst mode? (1=yes, 0=no)
% cLF            - apply Hipgh Pass Filter? ('none' or 'specify')
% HParam_in      - user defined cut off period for each session [1*nsess]
% cHF            - Low Pass Filter? ('none', 'Gaussian' or 'hrf')
% LParam_in      - user defined Gausssian FWHM (secs)
% cVi            - intrinsic correlations ('none' or 'AR(1)')
% bFcon          - Setup trial-specific F-contrasts? (1=yes, 0=no)


% =========================================================================
%function [xX,Sess] = spm_fmri_spm_ui
% Setting up the general linear model for fMRI time-series
% FORMAT [xX,Sess] = spm_fmri_spm_ui
%
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
% N.B.
% Burst Mode is a specialist design for intermittent epochs of acquisitions
% (used for example to allow for intercalated EEG recording).  Each burst
% is treated as a session but consistent within-session effects (e.g. T1
% effects) are modeled in X.bX.  The primary use of this mode is to generate
% parameter estimate images for a second level analysis.
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
% @(#)spm_fmri_spm_ui.m	2.26 Karl Friston, Jean-Baptiste Poline, Christian Buchel 99/05/19
%SCCSid  = '2.26';
SCCSid  = '';


%-GUI setup
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','fMRI stats model setup',0);
spm_help('!ContextHelp',mfilename)

global fmriDIR fmriTEST 

logfile = [fmriDIR filesep 'LOG' filesep 'analysis.log'];

str = ['cd ' fmriDIR filesep 'RESULTS/' kulRES]
eval(str);

%-Initialise output arguments in case return early
%xX   = [];
%Sess = [];


%	case 3
	% load pre-specified design matrix
	%---------------------------------------------------------------
	if sf_abort([2,3]), spm_clf(Finter), return, end
        str = ['load SPM_fMRIDesMtx'];
        eval(str);


	% get filenames
	%---------------------------------------------------------------
	nsess  = length(xX.iB);
	nscan  = zeros(1,nsess);
	for  i = 1:nsess
		nscan(i) = length(find(xX.X(:,xX.iB(i))));
	end
        P = [];

%	if nsess < 16
		for i = 1:nsess
			%str = sprintf('select scans for session %0.0f',i);
			%q   = spm_get(nscan(i),'.img',str);

                        q  = f99_P(DIRe(i,:),EXP(i,:));
 			P   = strvcat(P,q);
		end
%	else
%		str   = sprintf('select scans for this study');
%		P     = spm_get(sum(nscan),'.img',str);
%	end

	% Repeat time
	%---------------------------------------------------------------
	RT     = xX.RT;



% Assemble other deisgn parameters
%=======================================================================
spm_help('!ContextHelp',mfilename)
spm_input('Global intensity normalisation...',1,'d',mfilename)

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
%Global = spm_input(str,'+1','scale|none',{'Scaling' 'None'});
Global = {glob};

% Burst mode
%-----------------------------------------------------------------------
%if length(nscan) > 16 & ~any(diff(nscan))
%	spm_input('Burst mode?','+1','d',mfilename)
%	BM    = spm_input('Burst mode','+1','y/n',[1 0]);
%else
%	BM    = 0;
%end

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
%spm_input('Temporal autocorrelation options...','+1','d',mfilename)

% High-pass filtering
%-----------------------------------------------------------------------
%if BM
%	cLF = 'none';
%else
%	cLF = spm_input('High-pass filter?','+1','b','none|specify');
%end

switch cLF

        case 'specify'
	%---------------------------------------------------------------
%	HParam = spm_input(str,'+1','e',HParam,[1 nsess]);
        HParam = HParam_in;
        % LF description
	%---------------------------------------------------------------
	LFstr = sprintf('[min] Cutoff period %d seconds',min(HParam));

	case 'none'
	%---------------------------------------------------------------
	HParam = cell(1,nsess);
	LFstr  = cLF;
      
        case 'auto'
          cLF = 'specify';
        %---------------------------------------------------------------
	% default based on peristimulus time
	% param = cut-off period (max = 512, min = 32)
	%---------------------------------------------------------------
	HParam = 512*ones(1,nsess);
	for  i = 1:nsess
		for j = 1:length(Sess{i}.pst)
			HParam(i) = min([HParam(i) 2*max(RT + Sess{i}.pst{j})]);
		end
	end
	HParam = ceil(HParam);
        HParam_st = HParam;
	HParam(HParam < 32) = 32;
	str    = 'session cutoff period (secs)';
        % LF description
	%---------------------------------------------------------------
	LFstr = sprintf('[min] Cutoff period %d seconds',min(HParam));
        
        f99_log(logfile,['Default High Pass Filter is: ' mat2str(HParam)]);
end


% Low-pass filtering
%-----------------------------------------------------------------------
%if BM
%	cHF = 'none';
%else
%	cHF = spm_input('Low-pass filter?','+1','none|Gaussian|hrf');
%end

switch cHF

	case 'Gaussian'
	%---------------------------------------------------------------
%	LParam  = spm_input('Gaussian FWHM (secs)','+1','r',4);
        LParam  = LParam_in;
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
			'RT',		RT);
end


% intrinsic autocorrelations (Vi)
%-----------------------------------------------------------------------
str     = 'Model intrinsic correlations?';
cVimenu = {'none','AR(1)'};
%cVi     = spm_input(str,'+1','b',cVimenu);

%-Generate default trial-specific F-contrasts specified by session?
%-----------------------------------------------------------------------
%bFcon = spm_input('Setup trial-specific F-contrasts?','+1','y/n',[1,0],1,...
%		'batch',{},'trial_fcon');


% the interactive parts of spm_spm_ui are now finished: Cleanup GUI
%-----------------------------------------------------------------------
spm_clf(Finter);
spm('FigName','Configuring, please wait...',Finter,CmdLine);
spm('Pointer','Watch');


% Contruct K and Vi structs
%=======================================================================
K       = spm_filter('set',K);

% Adjust for missing scans
%-----------------------------------------------------------------------
%[xX,Sess,K,P,nscan,row] = spm_bch_tsampl(xX,Sess,K,P,nscan,row); %-SR

% create Vi struct
%-----------------------------------------------------------------------
Vi      = speye(sum(nscan));
xVi     = struct('Vi',Vi,'Form',cVi);
for   i = 1:nsess
	xVi.row{i} = row{i};
end


% get file identifiers and Global values
%=======================================================================
fprintf('%-40s: ','Mapping files') %-#
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

%-trial-specifc effects specified by Sess
%-----------------------------------------------------------------------
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
			'Interscan_interval',		sprintf('%0.2f',RT),...
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
spm_DesRep('DesMtx',xX,{VY.fname}',xsDes)
fprintf('%30s\n','...done')                                          %-#


%-Analysis Proper
%=======================================================================
spm_clf(Finter);
spm('FigName','fMRI stats models',Finter,CmdLine);
%if spm_input('estimate?',1,'b','now|later',[1,0],1)
calc = ['now'];
if strcmp(calc,'now')
	spm('Pointer','Watch')
   spm('FigName','Stats: estimating...',Finter,CmdLine);
   
   % Tor added this: load Vi
   disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
   disp(['Loading individual subject intrinsic autocorrelation'])
   disp(['Current Directory is ' pwd])
   disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
   % put Vi in file myVi in each subject directory.
   load myVi
   xX.xVi.Vi = Vi;
   
   
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
	abort = spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
	if abort, fprintf('%-40s: %30s\n\n',...
		'Abort...   (existing SPMstats files)',spm('time')), end
else
	abort = 0;
end
