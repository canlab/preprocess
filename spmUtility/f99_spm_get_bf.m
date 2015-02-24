function [BF,BFstr] = f99_spm_get_bf(name,T,dt,Fstr,Rov,Cov,pst_ev,h,m,HRF,TD,W,s,Sstr,rep,nsess,hrffunct)
% function [BF,BFstr] = spm_get_bf(name,T,dt,Fstr)
% creates basis functions for each trial type {i} in struct BF{i}
% FORMAT [BF BFstr] = spm_get_bf(name,T,dt,Fstr)
%
% name  - name{1 x n} name of trials or conditions
% T     - time bins per scan
% dt    - time bin length {seconds}
% Fstr  - Prompt string (usually indicates session)
%
% BF{i} - Array of basis functions for trial type {i}
% BFstr - description of basis functions specified
%_______________________________________________________________________
%
% spm_get_bf prompts for basis functions to model event or epoch-related
% responses.  The basis functions returned are unitary and orthonormal
% when defined as a function of peri-stimulus time in time-bins.
% It is at this point that the distinction between event and epoch-related 
% responses enters.
%_______________________________________________________________________
% @(#)spm_get_bf.m	2.10 Karl Friston 99/05/13

global fmriTEST

%-GUI setup
%-----------------------------------------------------------------------
spm_help('!ContextHelp',mfilename)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<4, Fstr = ''; end


% if no trials
%-----------------------------------------------------------------------
n      = length(name);
if ~n
	BF    = {};
	BFstr = 'none';
	return
end

% determine sort of basis functions
%-----------------------------------------------------------------------
Rtype = {'events',...
	 'epochs',...
	 'mixed'};
%if n == 1
%	Rtype = Rtype(1:2);
%	spm_input(name{1},1,'d',Fstr)
%else
%	spm_input(Fstr,1,'d')
%end
%Rov   = spm_input('are these trials',2,'b',Rtype);


for tr = 1:n        % loop over alle trial-types

    Rov_tr = Rov(tr,:);
  
    switch Rov_tr

	% assemble basis functions {bf}
	%===============================================================
	case 'ev'

	% model event-related responses
	%---------------------------------------------------------------
	Ctype = {
		'hrf (alone)',...
		'hrf (with time derivative)',...
		'hrf (with time and dispersion derivatives)',...
		'basis functions (Fourier set)',...
		'basis functions (Windowed Fourier set)',...
		'basis functions (Gamma functions)',...
		'basis functions (Gamma functions with derivatives)'};
	str   = 'Select basis set';
%	Cov   = spm_input(str,2,'m',Ctype);
        Cov_tr= Cov(tr);
	BFstr = Ctype{Cov_tr-4};


	% create basis functions
	%---------------------------------------------------------------
	if     Cov_tr == 8 | Cov_tr == 9

		% Windowed (Hanning) Fourier set
		%-------------------------------------------------------
		str   = 'window length {secs}';
%		pst   = spm_input(str,3,'e',32);
		pst_ev   = [0:dt:pst_ev]';
		pst_ev   = pst_ev/max(pst_ev);
%		h     = spm_input('order',4,'e',4);

		% hanning window
		%-------------------------------------------------------
		if Cov_tr == 8
			g = ones(size(pst_ev));
		else
			g = (1 - cos(2*pi*pst_ev))/2;
		end

		% zeroth and higher terms
		%-------------------------------------------------------
		bf    = g;
		for i = 1:h
			bf = [bf g.*sin(i*2*pi*pst_ev)];
			bf = [bf g.*cos(i*2*pi*pst_ev)];	
		end

	elseif Cov_tr == 10 | Cov_tr == 11


		% Gamma functions alone
		%-------------------------------------------------------
		pst_ev   = [0:dt:32]';
		dx    = 0.01;
		bf    = spm_gamma_bf(pst_ev);

		% Gamma functions and derivatives
		%-------------------------------------------------------
		if Cov_tr == 11
			bf  = [bf (spm_gamma_bf(pst_ev - dx) - bf)/dx];
		end


	elseif Cov_tr == 5 | Cov_tr == 6 | Cov_tr == 7


		% hrf and derivatives
		%-------------------------------------------------------
%                if ushrf
%                    if isempty(hrffunct)
%                    % 'with GUI'
%                       f_len = size(spm_hrf(dt));
%                       if spm_input('user specified hrf ?','+1','yes|no',[1 0])
%                          bf = spm_input(['user specified hrf function ' mat2str(f_len)],'+1');
%                       else
%                          [bf p] = spm_hrf(dt);
%                       end
%                    % 'without GUI'
%                    else
%                       bf = hrffunct{s};
%                    end
%                else
%                   [bf p] = spm_hrf(dt);
%                end

                [bf p] = spm_hrf(dt);
                if ~isempty(hrffunct)
                     bf = hrffunct(:);
                end

                % display of properties of the hrf function
                if fmriTEST & s==1 & tr==1
                     str = ['  The length of the hrf function is : ' ...
                           mat2str(length(spm_hrf(dt)))];
                     fprintf([str,' !!!\n\n']);
                end

                
		% add time derivative
		%-------------------------------------------------------
		if Cov_tr == 6 | Cov_tr == 7

			dp    = 1;
			p(6)  = p(6) + dp;
			D     = (bf(:,1) - spm_hrf(dt,p))/dp;
			bf    = [bf D(:)];
			p(6)  = p(6) - dp;

			% add dispersion derivative
			%-----------------------------------------------
			if Cov_tr == 7

				dp    = 0.01;
				p(3)  = p(3) + dp;
				D     = (bf(:,1) - spm_hrf(dt,p))/dp;
				bf    = [bf D(:)];
			end
		end
	end


	% Orthogonalize and fill in basis function structure
	%---------------------------------------------------------------
	bf    =  spm_orth(bf);
%	for i = 1:n
		BF{tr}  =  bf;
%	end


	% assemble basis functions {bf}
	%===============================================================
	case 'ep'


	% covariates of interest - Type
	%---------------------------------------------------------------
	Ctype = {'basis functions  (Discrete Cosine Set)',...
		 'basis functions  (Mean & exponential decay)',...
		 'fixed response   (Half-sine)',...
		 'fixed response   (Box-car)'};
	str   = 'Select type of response';
%	Cov   = spm_input(str,2,'m',Ctype);
        Cov_tr= Cov(tr);
        BFstr = Ctype{Cov_tr};


	% convolve with HRF?
	%---------------------------------------------------------------
	if Cov_tr == 1
		str = 'number of basis functions';
%		h   = spm_input(str,3,'e',2);
	end

	% convolve with HRF?
	%---------------------------------------------------------------
%	HRF   = spm_input('convolve with hrf',3,'b','yes|no',[1 0]);

	% ask for temporal differences
	%---------------------------------------------------------------
	str   = 'add temporal derivatives';
%	TD    = spm_input(str,4,'b','yes|no',[1 0]);
 

	% Assemble basis functions for each trial type
	%---------------------------------------------------------------
%        for i = 1:n

%		str   = ['epoch length {scans} for ' name{tr}];
%		W     = spm_input(str,'+1','r');
                %if s > 1
                if ~rep & nsess~=1
                   W_ses = W{s};
                else
                   W_ses = W{1};
                end
                W_tr  = W_ses(tr);
                %else
                %  W_tr  = W(tr);
                %end
		pst   = [1:W_tr*T]' - 1;
		pst   = pst/max(pst);

		% Discrete cosine set
		%-------------------------------------------------------
		if     Cov_tr == 1

			bf    = [];
			for j = 0:(m - 1)
				bf = [bf cos(j*pi*pst)];	
			end

		% Mean and exponential
		%-------------------------------------------------------
		elseif Cov_tr == 2
		
			bf    = [ones(size(pst)) exp(-pst/4)];

		% Half sine wave
		%-------------------------------------------------------
		elseif Cov_tr == 3

			bf    = sin(pi*pst);

		% Box car
		%-------------------------------------------------------
		elseif Cov_tr == 4

			bf    = ones(size(pst));

		end

		% convolve with hemodynamic response function - hrf
		%-------------------------------------------------------
		if HRF
			hrf   = spm_hrf(dt);
			[p q] = size(bf);
			D     = [];
			for j = 1:q
				D = [D conv(bf(:,j),hrf)];
			end
			bf    = D;
		end

		% add temporal differences if specified
		%-------------------------------------------------------
		if TD
			bf    = [bf [diff(bf); zeros(1,size(bf,2))]/dt];
		end

		% Orthogonalize and fill in Sess structure
		%-------------------------------------------------------
		BF{tr}         =  spm_orth(bf);

%	end


	% mixed event and epoch model
	%===============================================================
	case 'mixed'
	for i = 1:n

		BF(i)  = spm_get_bf(name(i),T,dt);

	end
	BFstr = 'mixed';

    end
end


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================


function bf = spm_orth(BF)
% recursive orthogonalization of basis functions
% FORMAT bf = spm_orth(bf)
%_______________________________________________________________________
bf    = BF(:,1);
bf    = bf/sqrt(sum(bf.^2));
for i = 2:size(BF,2)
	D     = BF(:,i);
	D     = D - bf*(pinv(bf)*D);
	if any(D)
		bf = [bf D/sqrt(sum(D.^2))];
	end
end


% compute Gamma functions functions
%-----------------------------------------------------------------------
function bf = spm_gamma_bf(u)
% returns basis functions used for Volterra expansion
% FORMAT bf = spm_gamma_bf(u);
% u   - times {seconds}
% bf  - basis functions (mixture of Gammas)
%_______________________________________________________________________
u     = u(:);
bf    = [];
for i = 2:4
        m   = 2^i;
        s   = sqrt(m);
        bf  = [bf spm_Gpdf(u,(m/s)^2,m/s^2)];
end
