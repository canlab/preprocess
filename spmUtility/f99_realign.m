function f99_realign(do_what,DIR,EXP,reslice_type,create_what,mask_images,adj_se,real_qual,rmflag)

% '==========================================================================='
%  ' REALIGNMENT AND RESLICING OF FUNCTIONAL IMAGES'
% '==========================================================================='
% '               --- The Script Parameters Explained ---                     '


%---------------------------------------------------------------------------
%   do_what? 1 = coregister only
%            2 = reslice only
%            3 = coregister and reslice
%   DIR      - [matrix van de directories v/d experimenten]
%   EXP      - [matrix van de experimentnamen]
%   reslice_type?
%   	      1  = Trilinear
%            -9  = Sinc
% 	     Inf = Fourrier (only if isotropic voxels)
%  create_what?
%            1 = All Images (1..n)| 2 = Images 2..n|'..
%	     3 = All Images + Mean Image| 4 = Mean Image Only'
%  mask_images? 
%            1 = yes, 0 = no
%  adj_se?   
%            1 = yes, 0 = no
%  registration quality? - number 1.00 to 0.001
%            'Quality 1.00  (slowest/most accurate) |Quality 0.90|' ...
%            'Quality 0.75|Quality 0.50|Quality 0.25|Quality 0.10|' ...
%            'Quality 0.05|Quality 0.01|' ...
%            'Quality 0.005|Quality 0.001 (fastest/poorest)'
%---OPTIONALY----------------------------------------------------------------
%   rmflags - status file in the LOG that controls if files can be
%             safely deleted
%---------------------------------------------------------------------------


% Within Mode Image Realignment
%___________________________________________________________________________
%
% This routine realigns a time-series of images acquired from the same subject 
% using a least squares approach and a 6 parameter (rigid body) spatial 
% transformation.  The first image in the list specified by the user is used
% as a reference to which all subsequent scans are realigned.  The reference
% scan does not have to the the first chronologically and it may be wise to
% chose a 'representative scan' in this role.
%
% For fMRI data an additional adjustment is made to the data that removes
% a tiny amount of the movement-related confounds of these effects.
% However, it may be preferable to include the functions of the estimated
% movement parameters as confounds in the statistics part.
%
%
% Uses
% Primarily to remove movement artefact in fMRI and PET time-series (or more 
% generally longitudinal studies)
%
%
% Inputs
% A series of *.img conforming to SPM data format (see 'Data Format').  The 
% relative displacement of the images should be small with respect to their 
% resolution.  This is usually easy to ensure for functional images (e.g. 
% fMRI, PET SPECT).
%
%
% Outputs
% The parameter estimation part writes out ".mat" files for each of the
% input images.  The part of the routine that writes the resliced images
% uses information in these ".mat" files and writes the realigned *.img
% files to the same subdirectory prefixed with an 'r' (i.e. r*.img).  The
% details of the transformation are displayed in the results window as
% plots of translation and rotation.
% A set of realignment parameters are saved for each session, named:
% realignment_params_*.txt.
%
%__________________________________________________________________________
% Refs:
%
% Friston KJ, Ashburner J, Frith CD, Poline J-B, Heather JD & Frackowiak
% RSJ (1995) Spatial registration and normalization of images Hum. Brain
% Map. 2:165-189
%
% Friston KJ, Williams SR, Howard R Frackowiak RSJ and Turner R (1995)
% Movement-related effect in fMRI time-series.  Mag. Res. Med. 35:346-355
%
% W. F. Eddy, M. Fitzgerald and D. C. Noll (1996) Improved Image
% Registration by Using Fourier Interpolation. Mag. Res. Med. 36(6):923-931
%
% R. W. Cox and A. Jesmanowicz (1999)  Real-Time 3D Image Registration
% for Functional MRI.  Submitted to MRM (April 1999) and avaliable from:
% http://varda.biophysics.mcw.edu/~cox/index.html.
%
%__________________________________________________________________________
%
% --- The Prompts Explained ---
%
% 'number of subjects'
% Enter the number of subjects you wish to realign.
%
% For fMRI, it will ask you the number of sessions for each subject.
% In the coregistration step, the sessions are first realigned to
% each other, by aligning the first scan from each session to the
% first scan of the first session.  Then the images within each session
% are aligned to the first image of the session.
% The parameter estimation is performed this way because it is assumed
% (rightly or not) that there may be systematic differences
% in the images between sessions.
% The adjustment step (correcting for resampling artifacts) is also
% performed completely independantly between each of the fMRI sessions.
%
% 'select scans for subject ..'
% Select the scans you wish to realign. All operations are relative
% to the first image selected.
%
% ......... Note that not all of the following prompts may be used: .........
%
% 'Which option?'
% 	'Coregister only'
% 	Only determine the parameters required to transform each of the
% 	images 2..n to the same space as image 1.
% 	The determined parameters for image XXXX.img are saved in the
%	file XXXX.mat. This is a Matlab file, containing the matrix 'M'.
% 	The location of an image voxel (in mm) can be determined by
% 	computing M(1:3,:)*[xcoord ycoord zcoord 1].
%	Note that if the coregistration is performed more than once on
%	the unresliced data, the starting estimates are obtained from
% 	parameters stored in the '.mat' files.
%	Note that for PET, the coregistration is a two step process.
%	First of all, the images are all realigned to the first in
%	the series.  A mean of these realigned images is created, and
%	a second pass realignment is performed to realign all the
%	images to the mean. Finally, the parameters are corrected
%	for any differences estimated by registering the first image in
%	the series to the mean image.
%
% 	'Reslice Only'
% 	Reslice the specified images according to the contents of the
% 	previously determined parameters. The images are resliced to be
% 	in the same space as the first one selected.  For fMRI, this is
%	the first image of the first session.
%
% 	'Coregister & Reslice'
% 	Combine the above two steps together.
%
%
% Options for reslicing:
%
% 'Create what?'
% 	'All Images (1..n)'
% 	This reslices all the images - including the first image selected
% 	- which will remain in it's original position.
%
%	'Images 2..n'
% 	Reslices images 2..n only. Useful for if you wish to reslice
% 	(for example) a PET image to fit a structural MRI, without
% 	creating a second identical MRI volume.
%
%	'All Images + Mean Image'
% 	In addition to reslicing the images, it also creates a mean of the
% 	resliced image.
%
%	'Mean Image Only'
% 	Creates the mean image only.
%
% 'Mask the images?'
% To avoid artifactual movement-related variance.
% Because of subject motion, different images are likely to have different
% patterns of zeros from where it was not possible to sample data.
% With masking enabled, the program searches through the whole time series
% looking for voxels which need to be sampled from outside the original
% images. Where this occurs, that voxel is set to zero for the whole set
% of images (unless the image format can represent NaN, in which case
% NaNs are used where possible).
%
% 'Adjust sampling errors?' (fMRI only)
% Adjust the data (fMRI) to remove interpolation errors arising from the
% reslicing of the data.  The adjustment for each fMRI session is performed
% independantly of any other session.  Bayesian statistics are used to
% attempt to regularize the adjustment in order to prevent an excessive
% amount of signal from being removed.  A priori variances for coefficients
% are assumed to be stationary and are estimated by translating the first
% image by a number of different distances using both Fourier and sinc
% interpolation.  This gives a ball park figure on how much error is
% likely to arise because of the approximations in the sinc interpolation.
% The certainty of the solution is obtained from the residuals after
% fitting the optimum linear combination of the basis functions through
% the data.  Estimates of certainty based on the residuals are
% unfortunately just an approximation.   
% We still don't fully understand the nature of the movement artifacts
% that arise using fMRI.  The current model is simply attempting to remove
% interpolation errors.  There are many other sources of error that the
% model does not attempt to remove.
% It is possible that adjusting the data without taking into account
% the design matrix for the statistics may be problematic when there are
% stimulous correlated movements, since adjusting seperately requires the
% assumption that the movements are independant from the paradigm.  It
% MAY BE BE BETTER TO INCLUDE THE ESTIMATED MOTION PARAMETERS AS CONFOUNDS
% WHEN THE STATISTICS ARE RUN.  The motion parameters are saved for each
% session, so this should be easily possible.
%
% 'Reslice Interpolation Method?'
% 	'Trilinear Interpolation'
% 	Use trilinear interpolation (first-order hold) to sample the images
%       during the writing of realigned images.
%
% 	'Sinc Interpolation'
% 	Use a sinc interpolation to sample the images during the writing
%	of realigned images.
% 	This is slower than bilinear interpolation, but produces better
% 	results. It is especially recommended for fMRI time series.
%	An 9x9x9 kernel is used to resample the images. 
%
%	'Fourier space Interpolation' (fMRI only)
%	Rigid body rotations are executed as a series of shears, which
%	are performed in Fourier space (Eddy et. al. 1996).  This routine
%	only supports cubic voxels (since zooms can not be done by
%	convolution in Fourier space).
%	No adjustment is available for this.
%
%__________________________________________________________________________
%
% The `.mat' files.
%
% This simply contains a 4x4 affine transformation matrix in a variable `M'.
% These files are normally generated by the `realignment' and
% `coregistration' modules.  What these matrixes contain is a mapping from
% the voxel coordinates (x0,y0,z0) (where the first voxel is at coordinate
% (1,1,1)), to coordinates in millimeters (x1,y1,z1).  By default, the
% the new coordinate system is derived from the `origin' and `vox' fields
% of the image header.
%  
% x1 = M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4)
% y1 = M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4)
% z1 = M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4)
%
% Assuming that image1 has a transformation matrix M1, and image2 has a
% transformation matrix M2, the mapping from image1 to image2 is: M2\M1
% (ie. from the coordinate system of image1 into millimeters, followed
% by a mapping from millimeters into the space of image2).
%
% These `.mat' files allow several realignment or coregistration steps to be
% combined into a single operation (without the necessity of resampling the
% images several times).  The `.mat' files are also used by the spatial
% normalisation module.
%__________________________________________________________________________
% @(#)spm_realign_ui.m	2.2 John Ashburner - with input from Oliver Josephs 00/01/17

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'realign',index_of_Analysis};
% or 
%    BCH.index0  = {'RealignCoreg',index_of_Analysis}; (when
%                   spm_realign is launched for edit_defaults 
%
%_______________________________________________________________________

global fmriTEST fmriDIR;

% fmri init stuff
%----------------
if nargin < 8 | nargin > 9
   error('This script needs 8 arguments (+ 1 optional rmflag)');
elseif nargin == 8
   rmflag = '';
end

D = [];
for i = 1:size(DIR,1)
   d = [deblank(fmriDIR) filesep deblank(DIR(i,:))];
   D = [D; d];
end
DIR = D;

if f99_exist(fmriDIR,'RESULTS') == 0
   str = ['!mkdir ' fmriDIR  filesep 'RESULTS'];
   eval(str);
end
if f99_exist(fmriDIR,'LOG') == 0
   str = ['!mkdir ' fmriDIR filesep 'LOG'];
   eval(str);
end
if f99_exist([fmriDIR filesep 'LOG'],'REALIGN.OK')
   str = ['!rm -f ' fmriDIR filesep 'LOG' filesep 'REALIGN.OK'];
   eval(str);
end

% Start logging
%--------------
logfile = [fmriDIR filesep 'LOG' filesep 'realign.log'];
tijd=spm('time');
disp('');
disp('*****************************************************************');
f99_log(logfile,['Realignment Script started at ' mat2str(tijd)]);
str = ['Start logging in ' logfile];
f99_log(logfile,str);

global MODALITY sptl_WhchPtn sptl_DjstFMRI sptl_CrtWht sptl_MskOptn SWD
global BCH sptl_RlgnQlty sptl_WghtRg;

%Setting the registration quality, which is a global value 
splt_RlgnQlty = real_qual;

%if (nargin == 0)
	% User interface.
	%_______________________________________________________________________
	%SPMid = spm('FnBanner',mfilename,'2.21');
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Realign');
	spm_help('!ContextHelp','spm_realign_ui.m');

	pos = 1;

	%n     = spm_input('number of subjects', pos, 'e', 1,...
        %                  'batch',{},'subject_nb');
	%if (n < 1)
	%	spm_figure('Clear','Interactive');
	%	return;
	%end
        n=1;
	P = cell(n,1);
	pos = pos + 1;
        MODALITY = 'FMRI';

        i=1;
	if strcmp(MODALITY,'FMRI'),
	    %ns = spm_input(['num sessions for subject 'num2str(i)], pos,...
            %    'e',1,'batch',{},'num_sessions');
            ns = size(DIR,1);
            pp = cell(1,ns);
            f99_log(logfile,' Coregistering: ');
            for s=1:ns,
                    p = '';
                    %while size(p,1)<1,
                    %	if isempty(BCH),
                    %		p = spm_get(Inf,'.img',...
                    %		['scans for subj ' num2str(i) ', sess' num2str(s)]);
                    %	else,
                    %		p = spm_input('batch',{'sessions',i},'images',s);
                    %	end;
                    %end;
                    p = f99_P(DIR(s,:),EXP(s,:));
                    nf  = size(p,1);
                    str = ['  SESSION ' spm_str_manip(DIR(s,:),'t') ' with number of timepoints = ' mat2str(nf)];
                    f99_log(logfile,str);
                    f99_checkP(DIR(s,:),EXP(s,:),1);
                    pp{s} = p;
            end
            P{i} = pp;

        else, %- no batch mode for 'PET'
            p  = cell(1,1);
            p{1} = '';
            %while size(p{1},1)<1,
            %      p{1} = spm_get(Inf,'.img',...
            %	  ['select scans for subject ' num2str(i)]);
            %end;
            p = f99_P(DIR(i,:),EXP(i,:));
            nf  = size(p,1);
            str = [' SESSION ' EXP(i,:) ' with number of timepoints = ' mat2str(nf)];
            f99_log(logfile,str);
            f99_checkP(DIR(i,:),EXP(i,:),1);
            P{i} = p;
        end;
%end;


	if strcmp(MODALITY,'PET'),
		FlagsC = struct('quality',sptl_RlgnQlty,'fwhm',8,'rtm',[]);
	else,
		FlagsC = struct('quality',sptl_RlgnQlty,'fwhm',6);
	end;

	if sptl_WhchPtn == 1,
		WhchPtn = 3;
	else,
%		WhchPtn = spm_input('Which option?', pos, 'm',...
%			'Coregister only|Reslice Only|Coregister & Reslice',...
%			[1 2 3],3,'batch',{},'option');
	        WhchPtn = do_what;
                pos = pos + 1;
	end;

	PW = '';
        sptl_WghtRg = 0; % NOT YET IMPLEMENTED !
	if (WhchPtn == 1 | WhchPtn == 3) & sptl_WghtRg,
		if spm_input(...
			['Weight the reference image(s)?'],...
			2, 'm',...
			['Dont weight registration|'...
			 'Weight registration'], [0 1], 1,...
			 'batch',{},'weight_reg'),

			if isempty(BCH),
				PW = spm_get(n,'.img',...
					'Weight images for each subj');
			else,
				PW = spm_input('batch',{'sessions',i},'weights',s);
			end;
		end;
	end;

	% Reslicing options
	%-----------------------------------------------------------------------
	if WhchPtn == 2 | WhchPtn == 3,
		FlagsR = struct('hold',1,'mask',0,'which',2,'mean',1);
		%FlagsR.hold = spm_input('Reslice interpolation method?',pos,'m',...
		%	     'Trilinear Interpolation|Sinc Interpolation|Fourier space Interpolation',...
		%	     [1 -9 Inf],2,'batch',{},'reslice_method');
		FlagsR.hold = reslice_type;
                strp = ['Trilinear Interpolation    ';...
                        'Sinc Interpolation         ';...
                        'Fourier space Interpolation'];
                if FlagsR.hold == 1
                    ff = 1;
                  elseif FlagsR.hold == -9
                    ff = 2;
                  elseif FlagsR.hold == Inf
                    ff = 3;
                end
                f99_log(logfile,[' Reslicing: ']);
                f99_log(logfile,['  using ',strp(ff,:)]);
                pos = pos + 1;

		%if sptl_CrtWht == 1,
		%	p = 3;
		%else
		%	p = spm_input('Create what?',pos,'m',...
		%		[' All Images (1..n)| Images 2..n|'...
		%		 ' All Images + Mean Image| Mean Image Only'],...
		%		[1 2 3 4],3,'batch',{},'create');
		%	pos = pos + 1;
		%end
		p = create_what;
                strp = ['All Images (1..n)       ';...
                        'Images 2..n             ';...
                        'All Images + Mean Image ';... 
                        'Mean Image Only         '];
                f99_log(logfile,['  ',deblank(strp(p,:)),' will be created']);      
   
                if (p == 1) FlagsR.which = 2; FlagsR.mean = 0; end
		if (p == 2) FlagsR.which = 1; FlagsR.mean = 0; end
		if (p == 3) FlagsR.which = 2; FlagsR.mean = 1; end
		if (p == 4) FlagsR.which = 0; FlagsR.mean = 1; end
		if FlagsR.which > 0,
                  
			%if sptl_MskOptn == 1,
			%	FlagsR.mask = 1;
			%else,
			%	if spm_input('Mask the resliced images?',pos,'y/n',...
                        %                      'batch',{},'mask') == 'y',
			%		FlagsR.mask = 1;
			%	end;
			%	pos = pos + 1;
			%end;
			if mask_images
                          FlagsR.mask = 1;
                          f99_log(logfile,['  Masking will be done']);
                        else
                          f99_log(logfile,['  Masking will not be done']);
                        end
      
                        
                        if strcmp(MODALITY, 'FMRI'),
				if finite(FlagsR.hold),
					%if sptl_DjstFMRI == 1,
					%	FlagsR.fudge = 1;
					%elseif sptl_DjstFMRI ~= 0,
					%	if spm_input(...
					%		'Adjust sampling errors?',pos,'y/n','batch',...
                                        %                {},'adjust_sampling_errors') == 'y',
					%	FlagsR.fudge = 1;
					%	end;
					%	pos = pos + 1;
					%end;
                                        if adj_se
                                          FlagsR.fudge = 1;
                                          f99_log(logfile,['  Adjust sampling errors will be done']);
                                        else
                                          f99_log(logfile,['  Adjust sampling errors will not be done']);
                                        end
				end
			end
		end
	end

        
if fmriTEST == 0
	spm('Pointer','Watch');
	for i = 1:n
		spm('FigName',['Realign: working on subject ' num2str(i)],Finter,CmdLine);
		fprintf('\rRealigning Subject %d: ', i);
		if WhchPtn==1 | WhchPtn==3,
			flagsC = FlagsC;
			if ~isempty(PW), flagsC.PW = deblank(PW(i,:)); end;
			spm_realign(P{i},flagsC);
		end
		if WhchPtn==2 | WhchPtn==3,
			spm_reslice(P{i},FlagsR)
		end;
	end
	fprintf('\r%60s%s', ' ',sprintf('\b')*ones(1,60));
	spm('FigName','Realign: done',Finter,CmdLine);
	spm('Pointer');
end

statusfile = [fmriDIR filesep 'LOG' filesep 'REALIGN.OK'];
f99_log(statusfile,'');
if ~isempty(rmflag)
   f99_rm(rmflag,DIR,EXP);
end

tijd=spm('time');
f99_log(logfile,['Realignment Script ended at ' mat2str(tijd)]);
disp('-----------------------------------------------------------------');
disp('');




