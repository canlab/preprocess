function f99_fmri_design_ma(RT,nscan,RES,STOC,rep,n_cond,cond_name,nulev,soa_stoc,oc_prob,prob_type,Sstr,soa_n_stoc,onset,Ptype,Pstr,Etype,h_exp,h_pol,sel,p_mod,Rov,Cov,pst_ev,h,m,HRF,TD,W,VOLT,regr_us,regr_us_na,durat,hrffunct)

% RT          - repitition time
% nscan       - number of scans per session
% RES         - directory where to write results
% STOC        - stochastic design (1=yes / 0=no)
% rep         - are sessions replicated exactly (1=yes / 0=no)
% n_cond           - number of conditions or trials
% cond_name   - matrix with names of conditions/trials
% 'stochastic design'
% s_STOC      - DO YOU WANT A STOCHASTIC DESIGN ('yes'/'no')
% ne          - null event (1=yes / 0=no)
% soa_stoc    - soa for each trial type
% oc_prob     - matrix with occurence probabilities for each event
% prob_type   - stationary (1) or modulated (0) occurence probabilities
% 'non stochastic design'
% N_STOC      - DO YOU WANT A DETERMINISTIC DESIGN ('yes'/'no')
% Sstr        - type of SOA ('fixed' or 'variable')
% onset       - time to first trial (scans)
% durat       - variable durations for all sessions and conditions 
%                   (ncond x # onsets for each cond) for each session 
% soa_n_stoc  - soa for each trial type (only if 'fixed')
% 'modulation'
% Ptype       - parametric modulation ('none','time','other')
% Pstr        - name of parameter (in case of 'other')
% Etype       - expansion type ('linear','exponen','polynom')
% h_exp       - time constant (in case of 'exponen')
% h_pol       - order of polynomial expansion (in case of 'polynom')
% sel         - 
% p_mod       - matrix with parameters for each selected covariate
% 'basisfunctions to use'
% Rov         - matrix of types of trials ('ev' or 'ep')
%                * event ('ev') or epoch ('ep')
%                * one string per trial type (e.g. ['ep';'ep';'ev'])
% hrffunct    - user specified hrf function of this subject
%                length of this function can be checked using fmriTEST
% Cov         - model type for covariates of interest
%                1 = basis functions  (Discrete Cosine Set
%                2 = basis functions  (Mean & exponential decay)
%                3 = fixed response   (Half-sine)
%                4 = fixed response   (Box-car)
%             - model type for event-related responses
%                5 = hrf (alone)
%                6 = hrf (with time derivative)
%                7 = hrf (with time and dispersion derivatives)
%                8 = basis functions (Fourier set)
%                9 = basis functions (Windowed Fourier set)
%               10 = basis functions (Gamma functions)
%               11 = basis functions (Gamma functions with derivatives)
% pst_ev      - window length (for Cov=8 Cov=9)        'event'
% h           - order of Fourier set (for Cov=8 Cov=9)
% m           - number of basis functions (if Cov=1)   'epoch'
% HRF         - convolve with HRF (1=yes / 0=no)
% TD          - add temporal differences (1=yes / 0=no)
% W           - epoch length (scans) for each condition  
% VOLT        - interactions among trials (Volterra series)? (1=yes / 0=no)
% regr_us     - matrix with user specified regressors
% regr_us_na  - matrix with strings of names of regr_us


% ======================================================================
% function [xX,Sess] = spm_fMRI_design(nscan,RT)
% Assembles a design matrix for fMRI studies
% FORMAT [xX,Sess] = spm_fMRI_design(nscan,RT)
%
% nscan   - n vector {nscan(n) = number of scans in session n}
% RT      - intercans interval {seconds}
%
% xX            - structure describing design matrix
% xX.X          - design matrix
% xX.dt         - time bin {secs}
% xX.RT         - Repetition time {secs}
% xX.iH         - vector of H partition (condition effects)      indices,
% xX.iC         - vector of C partition (covariates of interest) indices
% xX.iB         - vector of B partition (block effects)          indices
% xX.iG         - vector of G partition (nuisance variables)     indices
% xX.Xnames     - cellstr of effect names corresponding to columns
%                 of the design matrix
%
% Sess{s}.BFstr    - basis function description string
% Sess{s}.DSstr    - Design description string
% Sess{s}.rep      - session replication flag
% Sess{s}.row      - scan   indices      for session s
% Sess{s}.col      - effect indices      for session s
% Sess{s}.name{i}  - of ith trial type   for session s
% Sess{s}.ind{i}   - column indices      for ith trial type {within session}
% Sess{s}.bf{i}    - basis functions     for ith trial type
% Sess{s}.sf{i}    - stick functions     for ith trial type
% Sess{s}.ons{i}   - stimuli onset times for ith trial type (secs)
% Sess{s}.pst{i}   - peristimulus times  for ith trial type (secs)
% Sess{s}.Pv{i}    - vector of paramters for ith trial type
% Sess{s}.Pname{i} - name   of paramters for ith trial type
%
% saves SPM_fMRIDesMtx.mat (xX Sess)
%_______________________________________________________________________
%
% spm_fMRI_design allows you to build design matrices with separable
% session-specific partitions.  Each partition may be the same (in which
% case it is only necessary to specify it once) or different.  Responses
% can be either event- or epoch related, where the latter model prolonged
% and possibly time-varying responses to state-related changes in
% experimental conditions.  Event-related response are modelled in terms
% of responses to instantaneous events.  Mathematically they are both
% modelled by convolving a series of delta (or stick) functions,
% indicating the onset of an event or epoch with a set of basis
% functions.  These basis functions can be very simple, like a box car,
% or may model voxel-specific forms of evoked responses with a linear
% combination of several basis functions (e.g.  a Fourier set).  Basis
% functions can be used to plot estimated responses to single events or
% epochs once the parameters (i.e.  basis function coefficients) have
% been estimated.  The importance of basis functions is that they provide
% a graceful transition between simple fixed response models (like the
% box-car) and finite impulse response (FIR) models, where there is one
% basis function for each scan following an event or epoch onset.  The
% nice thing about basis functions, compared to FIR models, is that data
% sampling and stimulus presentation does not have to be sychronized
% thereby allowing a uniform and unbiased sampling of peri-stimulus time.
%
% spm_fMRI_design allows you to combine both event- and epoch-related
% responses in the same model.  You are asked to specify the number
% of condition (epoch) or trial (event) types.  Epoch and event-related 
% responses are modeled in exactly the same way by first specifying their
% onsets [in terms of stimulus onset asynchronies (SOAs) or explicit onset 
% times] and then convolving with appropriate basis functions (short ones
% for event-related models and longer ones for epoch-related respones).
% Enter 0 to skip these if you only want to use regressors you have designed 
% outside spm_fMRI_design).
%
% Interactions or response modulations can enter at two levels.  Firstly
% the stick function itself can be modulated by some parametric variate
% (this can be time or some trial-specific variate like reaction time)
% modeling the interaction between the trial and the variate or, secondly
% interactions among the trials themselves can be modeled using a Volterra
% series formulation that accommodates interactions over time (and therefore
% within and between trial types).  The first sort of interaction is
% specified by extra (modulated) stick functions in Sess{s}.sf{i}.  If
% a polynomial expansion of the specified variate is requested there will
% be more than one additional column.  The corresponding name of the
% explanatory variables in X.Xname is Sn(s) trial i(p)[q] for the qth
% order expansion of the variate convolved with the pth basis function
% of the ith trial in the sth session.  If no parametric variate is
% specified the name is simply Sn(s) trial i(p).  Interactions among
% and within trials enter as new trial types but do not have .pst or .ons
% fields.  These interactions can be characterized later, in results, in
% terms of the corresponding second order Volterra Kernels.
%
% The design matrix is assembled on a much finer time scale (X.dt) than the 
% TR and is then subsampled at the acquisition times.
%
% Sess{s}.ons{i} contains stimulus onset times in seconds relative to the
% timing of the first scan and are provided to contruct stimuli when using 
% stochastic designs
%
%
% Notes on spm_get_ons, spm_get_bf and spm_Volterra are included below
% for convenience.
%
%                           ----------------
%
% spm_get_ons contructs a cell of sparse delta functions specifying the
% onset of events or epochs (or both). These are convolved with a basis set
% at a later stage to give regressors that enter into the design matrix.
% Interactions of evoked responses with some parameter (time or a specified 
% variate Pv) enter at this stage as additional columns in sf with each delta
% function multiplied by the [expansion of the] trial-specific parameter.
% If parametric modulation is modeled, P contains the original variate and
% Pname is its name.  Otherwise P{i} = [] and Pname{i} = '';
%
%                           ----------------
%
% spm_get_bf prompts for basis functions to model event or epoch-related
% responses.  The basis functions returned are unitary and orthonormal
% when defined as a function of peri-stimulus time in time-bins.
% It is at this point that the distinction between event and epoch-related 
% responses enters.
%
%                           ----------------
%
% For first order expansions spm_Volterra simply convolves the causes
% (e.g. stick functions) in SF by the basis functions in BF to create
% a design matrix X.  For second order expansions new entries appear
% in IND, BF and name that correspond to the interaction among the
% orginal causes (if the events are sufficiently close in time).
% The basis functions for these are two dimensional and are used to
% assemble the second order kernel in spm_graph.m.  Second order effects
% are computed for only the first column of SF.
%
%_______________________________________________________________________
% @(#)spm_fMRI_design.m	2.17 Karl Friston 99/05/19
SCCSid  = '2.17';

% fmri init stuff
%----------------
global fmriDIR fmriTEST functie

if f99_exist(fmriDIR,'RESULTS') == 0
  str = ['!mkdir ' fmriDIR '/RESULTS'];
  eval(str);
end

% Change to results directory
%----------------------------
str = [fmriDIR '/RESULTS/'];
if f99_exist(str,RES) == 0
  str = ['!mkdir ' str RES];
  eval (str);
end

str = ['cd ' fmriDIR '/RESULTS/' RES];
str = deblank(str);
eval (str)


%-GUI setup
%-----------------------------------------------------------------------
SPMid = spm('SFnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','fMRI stats model setup',0);
spm_help('!ContextHelp',mfilename)


% construct Design matrix {X} - cycle over sessions
%=======================================================================

% Initialize variables
%-----------------------------------------------------------------------
%STOC   = 0;

% global parameters
%-----------------------------------------------------------------------
global fMRI_T; 
global fMRI_T0; 

if isempty(fMRI_T),  fMRI_T  = 16; end;
if isempty(fMRI_T0), fMRI_T0 = 16; end;

% get nscan and RT if no arguments
%-----------------------------------------------------------------------
%if nargin < 2
%	spm_input('Basic parameters...',1,'d',mfilename)
%	RT     = spm_input('Interscan interval {secs}','+1','r',[],1);
%end
%if nargin < 1
%	nscan  = spm_input(['scans per session e.g. 64 64 64'],'+1');
%	STOC   = 1;
%end
nsess  = length(nscan);
dt     = RT/fMRI_T;					% time bin {secs}


% separate specifications for non-relicated sessions
%-----------------------------------------------------------------------
%rep = 0;
%if nsess > 1 & ~any(nscan - nscan(1))
%	rep = spm_input(['are sessions replicated exactly'],'+1','yes|no',[1 0]);
%end
Xx    = [];
Xb    = [];
Xname = {};
Bname = {};
Sess  = {};
for s = 1:nsess

	if rep
		Fstr='All sessions';
	else
		Fstr = sprintf('Session %d/%d',s,nsess);
	end


	% specify design for this session?
	%---------------------------------------------------------------
        if   (s == 1) | ~rep
	X       = [];				% design matrix
	B       = [];				% block partition
	D       = [];				% covariates
	PST     = {};				% Peri-stimulus times
	ONS     = {};				% onset times
	IND     = {};				% design matrix indices
	name    = {};				% condition names
	Xn      = {};				% regressor names

	% Event/epoch related responses			
	%===============================================================
	k       = nscan(s);

	% event/epoch onsets {ons} and window lengths {W}
	%---------------------------------------------------------------
   
   [SF,name,Pv,Pname,DSstr] = f99_spm_get_ons(k,fMRI_T,dt,STOC,Fstr,n_cond,cond_name,nulev,soa_stoc,oc_prob,prob_type,Sstr,soa_n_stoc,onset,Ptype,Pstr,Etype,h_exp,h_pol,sel,p_mod,s,rep,nsess,durat);

        ntrial     = size(SF,2);


	% get basis functions for responses
	%---------------------------------------------------------------
	[BF BFstr] = f99_spm_get_bf(name,fMRI_T,dt,Fstr,Rov,Cov,pst_ev,h,m,HRF,TD,W,s,Sstr,rep,nsess,hrffunct);


	%-Reset ContextHelp
	%---------------------------------------------------------------
	spm_help('!ContextHelp',mfilename)


	% peri-stimulus {PST} and onset {ONS} times in seconds
	%---------------------------------------------------------------
   	   	for   i = 1:ntrial
       on     = find(SF{i}(:,1))*dt;
		pst    = [1:k]*RT - on(1);
		for  j = 1:length(on)
			u      = [1:k]*RT - on(j);
			v      = find(u >= -1);
			pst(v) = u(v);
		end
		ONS{i} = on;
		PST{i} = pst;
	end

%	spm_input('Additional design matrix options...',1,'d',mfilename)

	% convolve with basis functions
	%---------------------------------------------------------------

        if ntrial
%		str    = 'interactions among trials (Volterra)';
%		if spm_input(str,'+1','y/n',[1 0]);
                if VOLT
			[X Xn IND BF name] = spm_Volterra(SF,BF,name,2);
		else
			[X Xn IND]         = spm_Volterra(SF,BF,name,1);
		end

		% Resample design matrix {X} at acquisition times
		%-------------------------------------------------------
		X     = X([0:k-1]*fMRI_T + fMRI_T0,:);
	end




	% get user specified regressors
	%===============================================================
%	c     = spm_input('user specified regressors','+1','w1',0);
        if size(regr_us)
            if ~rep & nsess~=1
                ses_regr_us = regr_us{s};
            else
                ses_regr_us = regr_us{1};
            end
            c = size(ses_regr_us);
            while size(D,2) < c
                str   = sprintf('regressor %i',size(D,2) + 1);
                %		D     = [D spm_input(str,2,'e',[],[k Inf])];
                D     = [D ses_regr_us(size(D,2)+1,:)'];
            end
            if      c & length(DSstr)
                DSstr = [DSstr '& user specified covariates '];
            elseif  c
                DSstr = 'User specified covariates ';
            end

	% append regressors and names
        %---------------------------------------------------------------
            if ~rep & nsess~=1
                ses_regr_us_na = regr_us_na{s};
            else
                ses_regr_us_na = regr_us_na{1};
            end
            for i = 1:size(D,2)
                X           = [X D(:,i)];
                str         = sprintf('regressor: %i',i);
                Xn{end + 1} = ses_regr_us_na(i,:);
            end
        end
            

	% Confounds: Session effects 
	%===============================================================
	B      = ones(k,1);
	Bn{1}  = sprintf('constant');

	% end specification
	%---------------------------------------------------------------
	end

	% Session structure
	%---------------------------------------------------------------
	Sess{s}.BFstr = BFstr;
	Sess{s}.DSstr = DSstr;
	Sess{s}.rep   = rep;
	Sess{s}.row   = size(Xx,1) + [1:k];
	Sess{s}.col   = size(Xx,2) + [1:size(X,2)];
	Sess{s}.name  = name;
	Sess{s}.ind   = IND;
	Sess{s}.bf    = BF;
	Sess{s}.sf    = SF;
	Sess{s}.pst   = PST;
	Sess{s}.ons   = ONS;
	Sess{s}.Pv    = Pv;
	Sess{s}.Pname = Pname;

	% Append names
	%---------------------------------------------------------------
	q     = length(Xname);
	for i = 1:length(Xn) 
		Xname{q + i}  = [sprintf('Sn(%i) ',s) Xn{i}];
	end
	q     = length(Bname);
	for i = 1:length(Bn) 
		Bname{q + i}  = [sprintf('Sn(%i) ',s) Bn{i}];
	end

	% append into Xx and Xb
	%===============================================================
	[x y]   = size(Xx);
	[i j]   = size(X);
	Xx(([1:i] + x),([1:j] + y)) = spm_detrend(X);
	[x y]   = size(Xb);
	[i j]   = size(B);
	Xb(([1:i] + x),([1:j] + y)) = B;
end

% finished
%-----------------------------------------------------------------------
xX     = struct(	'X',		[Xx Xb],...
			'dt',		dt,...
			'RT',		RT,...
			'iH',		[],...
			'iC',		[1:size(Xx,2)],...
			'iB',		[1:size(Xb,2)] + size(Xx,2),...
			'iG',		[],...
			'Xnames',	{[Xname Bname]});


%-End: Save SPM_fMRIDesMtx.mat
%-----------------------------------------------------------------------
fprintf('\t%-32s: ','Saving fMRI design')                            %-#
save SPM_fMRIDesMtx SPMid xX Sess
fprintf('%30s\n\n','...SPM_fMRIDesMtx.mat saved')                    %-#
spm_input('!DeleteInputObj')


% 'review a specified model'
%---------------------------------------------------------------
spm_clf(Finter)
[xX,Sess] = spm_fMRI_design_show(xX,Sess);
return
