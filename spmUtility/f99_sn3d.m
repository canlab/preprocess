function f99_sn3d(option,objima,mask_object,maskima,DIR,EXP,templima,mask_templ,spec_brainmask,...
   affinet,ase,nbf,nli,nlr,intplt,bbt,voxt,rmflag)

% '==========================================================================='
%  ' SPATIAL (STEREOTATIC) NORMALIZATION'
% '==========================================================================='
% '               --- The Script Parameters Explained ---                     '


%---------------------------------------------------------------------------
%  option - 'Do What?'
% 	 1 = 'Determine Parameters Only'
%        2 = 'Write Normalised Only'
% 	 3 = 'Determine Parameters & Write Normalised'
%  objima - 'Image to determine parameters from (objectimage), OR
%        in case of write normalised only, image for which parameters exist
%        e.g. 'anat/anatomy.img'
%  mask_ojb - 'Mask the object?'
%        1 = yes, 0 = no
%  maskima  - 'image that contains the mask for the objectimage'
%        e.g. 'anat/noskull.img'
%  DIR - '[matrix of dirs]'
%  EXP - '[matrix of experiments]'
%  templima - '[one or more templates]'
%       The image will be fitted (in the least 
%       squares sense) to the optimum linear combination of these templates.
%  mask_templ 
%       0 = 'no masking'
%       1 = 'mask the template image with the mask in SPM/apriori directory'
%       2 = 'mask with the spec_brainmask (user defined)'
%  spec_brainmask - 'user specified brainmask for template'
%  affinet - 'Affine Normalisation Type?'
%       'rigid1' - rigid body normalisation
%       'rigid2'
%       'rigid3'
%       'affine1'
%       'affine2'
%       'affine3'  - this is default SPM99
%  ase - 'Affine Starting Estimates?'
%	'Neurological Convention (R is R)'
% 	The starting estimates do not include a left right flip.
% 	[0 0 0 0 0 0  1 1 1 0 0 0].
%	'Radiological Convention (L is R)'
% 	[0 0 0 0 0 0 -1 1 1 0 0 0].
%  nbf - '# Basis Functions (x y z)'
%       1=none  -  USE THIS IF AFFINE ONLY NORMALISATION WANTED
%       2=2x2x2 | 3=2x3x2 | 4=3x3x3 | 5=3x4x3 | 6=4x4x4 | 7=4x5x4 | 8=5x5x5
%       9=5x6x5 | 10=6x6x6 | 11=6x7x6| 12=7x7x7 | 13=7x8x7 | 14=8x8x8
%  nli - '# Nonlinear Iterations?'
%       can be set to : 1, 3, 5, 8, 12 or 16
%  nlr - 'Nonlinear Regularization'
%       can be set to : 1, 0.1, 0.01, 0.001 or 0.0001
%	from 'Extremely Heavy regularization' to 'Very Light regularization'
%       0.01 is SPM99 default, use Heavy is distorted
%  intplt - 'Interpolation Method?'
% 	0=Nearest Neighbour | 1=Bilinear Interpolation |...
%	-9=Sinc Interpolation (9x9x9)
%  bbt - 'Bounding Box?'
%       [-x x -y y -z z] in mm e.g. [-78 78 -112 76  -50 85]  (Default SPM99)
%  voxt - 'Voxel Sizes '
%       The voxel sizes (x, y & z, in mm) of the written normalised images.
%       [x y z] in mm e.g. [2 2 2]
%---OPTIONALY----------------------------------------------------------------
%   rmflags - status file in the LOG that controls if files can be
%             safely deleted
%---------------------------------------------------------------------------


%function spm_sn3d(P,matname,bb,Vox,params,spms,brainmask,objmask)
% Spatial (stereotactic) normalization
% ___________________________________________________________________________
% 
% This module spatially (stereotactically) normalizes MRI, PET or SPECT
% images into a standard space defined by some ideal model or template
% image[s].  The template images supplied with SPM conform to the space
% defined by the ICBM, NIH P-20 project, and approximate that of the
% the space described in the atlas of Talairach and Tournoux (1988).
% The transformation can also be applied to any other image that has
% been coregistered with these scans.
% 
% 
% Mechanism
% Generally, the algorithms work by minimising the sum of squares
% difference between the image which is to be normalised, and a linear
% combination of one or more template images.  For the least squares
% registration to produce an unbiased estimate of the spatial
% transformation, the image contrast in the templates (or linear
% combination of templates) should be similar to that of the image from
% which the spatial normalization is derived.  The registration simply
% searches for an optimum solution.  If the starting estimates are not
% good, then the optimum it finds may not find the global optimum.
% 
% The first step of the normalization is to determine the optimum
% 12-parameter affine transformation.  Initially, the registration is
% performed by matching the whole of the head (including the scalp) to
% the template.  Following this, the registration proceeded by only
% matching the brains together, by appropriate weighting of the template
% voxels.  This is a completely automated procedure (that does not
% require ``scalp editing'') that discounts the confounding effects of
% skull and scalp differences.   A Bayesian framework is used, such that
% the registration searches for the solution that maximizes the a
% posteriori probability of it being correct.  i.e., it maximizes the
% product of the likelihood function (derived from the residual squared
% difference) and the prior function (which is based on the probability
% of obtaining a particular set of zooms and shears).
% 
% The affine registration is followed by estimating nonlinear deformations,
% whereby the deformations are defined by a linear combination of three
% dimensional discrete cosine transform (DCT) basis functions.  The default
% options result in each of the deformation fields being described by 1176
% parameters, where these represent the coefficients of the deformations in
% three orthogonal directions.  The matching involved simultaneously
% minimizing the membrane energies of the deformation fields and the
% residual squared difference between the images and template(s).
% 
% An option is provided for allowing weighting images (consisting of pixel
% values between the range of zero to one) to be used for registering
% abnormal or lesioned brains.  These images should match the dimensions
% of the image from which the parameters are estimated, and should contain
% zeros corresponding to regions of abnormal tissue.
% 
% 
% Uses
% Primarily for stereotactic normalization to facilitate inter-subject
% averaging and precise characterization of functional anatomy.  It is
% not necessary to spatially normalise the data (this is only a
% pre-requisite  for  intersubject averaging or reporting in the
% Talairach space).  If you wish to circumnavigate this step  (e.g. if
% you have single slice data or do not have an appropriate high
% resolution MRI scan) simply specify where you think the  anterior
% commissure  is  with  the  ORIGIN in the header of the first scan
% (using the 'Display' facility) and proceed directly  to  'Smoothing'
% or 'Statistics'.
% 
% 
% Inputs
% The first input is the image which is to be normalised. This image
% should be of the same modality (and MRI sequence etc) as the template
% which is specified. The same spatial transformation can then be
% applied to any other images of the same subject.  These files should
% conform to the SPM data format (See 'Data Format'). Many subjects can
% be entered at once, and there is no restriction on image dimensions
% or voxel size.
% 
% Providing that the images have a correct ".mat" file associated with
% them, which describes the spatial relationship between them, it is
% possible to spatially normalise the images without having first
% resliced them all into the same space. The ".mat" files are generated
% by "spm_realign" or "spm_coregister".
% 
% Default values of parameters pertaining to the extent and sampling of
% the standard space can be changed, including the model or template
% image[s].
% 
% 
% Outputs
% All normalized *.img scans are written to the same subdirectory as
% the original *.img, prefixed with a 'n' (i.e. n*.img).  The details
% of the transformations are displayed in the results window, and the
% parameters are saved in the "*_sn3d.mat" file.
% 
% 
%____________________________________________________________________________
% Refs:
% K.J. Friston, J. Ashburner, C.D. Frith, J.-B. Poline,
% J.D. Heather, and R.S.J. Frackowiak
% Spatial Registration and Normalization of Images.
% Human Brain Mapping 2:165-189(1995)
% 
% J. Ashburner, P. Neelin, D.L. Collins, A.C. Evans and K. J. Friston
% Incorporating Prior Knowledge into Image Registration.
% NeuroImage 6:344-352 (1997)
%
% J. Ashburner and K. J. Friston
% Nonlinear Spatial Normalization using Basis Functions.
% Human Brain Mapping 7(4):in press (1999)
%
%_______________________________________________________________________
%
% --- The Prompts Explained ---
%
% 'Which option?'
% 	'Determine Parameters Only'
% 	Computes the parameters which best fit the image to the
% 	template(s),
% 	and saves them to the file imagename'_sn3d.mat'.
%
% 	'Write Normalised Only'
% 	Select the appropriate '_sn3d.mat' file, and the images you
% 	wish to have normalised. The normalised images will be written
% 	prefixed with 'r'.
%
% 	'Determine Parameters & Write Normalised'
% 	Combines the above two steps into one.
%
% Options for Determine Parameters:
%
% 'select Template(s) '
% Select one or more templates. The image will be fitted (in the least 
% squares sense) to the optimum linear combination of these templates.
%
% 'Normalisation Type?'
% 	'Default Normalisation'
% 	This is a 12 parameter affine normalisation, followed by
% 	5 iterations of the nonlinear normalisation using 4x5x4 basis
% 	functions. If this doesn't work, or if you wish to push the
% 	normalisation a bit harder, try a custom normalisation.
% 	The default arguments for the custom normalisation have
% 	(D) next to them.
%
% 	'Affine Only'
% 	Only perform an affine ('brain in a box') normalisation.
% 	The affine normalisation corrects for gross differences in
% 	brain size, and position.
%
% 	'Custom Affine & Nonlinear'
% 	Perform an affine and also a 3D nonlinear normalisation. The
% 	nonlinear part of the normalisation corrects for more subtle
% 	differences in brain shape.
%
% 'Affine Starting Estimates?'
% The normalised images should use neurological convention for their
% orientation. If any flipping needs doing, it should be done at this
% stage.
%
%	'Neurological Convention (R is R)'
% 	The starting estimates do not include a left right flip.
% 	[0 0 0 0 0 0  1 1 1 0 0 0].
%
%	'Radiological Convention (L is R)'
% 	The starting estimates flip left right (since the template uses
% 	Neurological Convention).
% 	[0 0 0 0 0 0 -1 1 1 0 0 0].
%
% 	'Custom Starting Estimates'
% 	Enter starting estimates in the following order:
% 		x translation (mm)
% 		y translation (mm)
% 		z translation (mm)
% 		x rotation about - {pitch} (radians)
% 		y rotation about - {roll}  (radians)
% 		z rotation about - {yaw}   (radians)
% 		x scaling
% 		y scaling
% 		z scaling
% 		x affine
% 		y affine
% 		z affine
%	For images aquired in different orientations, it is possible
% 	to use starting estimates which reflect these orientations.
% 	eg. a pitch of +/- pi/2 radians for coronal images.
% 	    a roll  of +/- pi/2 radians for saggital images.
% 	For volumes which are flipped, then the appropriate scaling
% 	can be set to -1 in the starting estimates.
%
% '# Nonlinear Iterations?'
% Try about 12 iterations.
%
% '# Basis Functions (x y z)'
% The deformation field is computed from a linear combination of (3D)
% basis images. The basis images used are those which make up the
% lowest frequency components of a 3D discrete cosine transform.
% What is entered here is the dimensions of this transform.
%
% 'Nonlinear Regularization'
% Pick a medium value.  However, if your normalized images appear
% distorted, then it may be an idea to increase the amount of
% regularization - or even just use an affine normalization.
% The regularization influences the smoothness of the deformation
% fields.
%
% 'Mask brain when registering?'
% Applies a weighting mask to the template(s) during the parameter
% estimation.  By default, weights in and around the brain have
% values of one whereas those clearly outside the brain have values
% of zero.  This is an attempt to base the normalization purely upon
% the shape of the brain, rather than the shape of the head (since
% low frequency basis functions can not really cope with variations
% in skull thickness).
% The option is now available for a user specified weighting image.
% This should have the same dimensions and mat file as the template
% images, with values in the range of zero to one.
%
% 'Mask object brain when registering?'
% Applies a weighting mask to the object image(s) during the parameter
% estimation.  Weights are as for the template mask (0-1).  Used
% (usually) to prevent unusual or abnormal areas of brain (e.g.
% stroke, tumour) influencing normalisation to normal brain
%
% Options for Write Normalised:
% 'select Normalisation Parameter Set'
% Select the '_sn3d.mat' file.
%
% 'Interpolation Method?'
% The method by which the images are sampled when being written in a
% different space.
% 	'Nearest Neighbour'
% 		- Fastest, but not normally recommended.
% 	'Bilinear Interpolation'
% 		- OK for PET, or realigned fMRI.
%
% 	'Sinc Interpolation'
%		- With sinc interpolation, a sinc function with a
%		  Hanning filter envelope is used to resample the data
%		  (9x9x9 kernel).
%
% 'Bounding Box?'
% The bounding box (in mm) of the volume which is to be written
% (relative to the anterior commissure).
%
% 'Voxel Sizes '
% The voxel sizes (x, y & z, in mm) of the written normalised images.
%
%_______________________________________________________________________
% @(#)spm_sn3d.m	2.26	John Ashburner MRCCU/FIL 99/10/29
% With suggested modifications by Matthew Brett of the MRCCU

% Programmers notes
%-----------------------------------------------------------------------
%
% FORMAT spm_sn3d('Defaults')
% acts as a user interface for setting normalisation defaults.
%
%-----------------------------------------------------------------------
%
% FORMAT spm_sn3d(P,matname,bb,Vox,params,spms,brainmask, objmask)
% P         - image(s) to normalize
% matname   - name of file to store deformation definitions
% bb        - bounding box for normalized image
% Vox       - voxel sizes for normalized image
% params(1) - number of basis functions in X
% params(2) - "      "  "     "         "  Y
% params(3) - "      "  "     "         "  Z
% params(4) - number of iterations for elastic deformation
%             Setting any of these parameters to 0 will force
%             the program to perform only the affine
%             normalization.
% params(5) - smoothing for image (mm).
% params(6) - amount of regularization.  Higher values produce smoother deformation
%             fields.
% spms      - template image(s).
% brainmask - Weighting image for template(s)
% objmask   - Weighting image for object images



% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'normalize',index_of_Analysis};
% or 
%    BCH.index0  = {'normalisation',index_of_Analysis}; (when
%                   spm_spn3d is launched for edit_defaults 
%
%_______________________________________________________________________

global SWD sptl_Vx sptl_BB sptl_NBss sptl_Ornt sptl_CO sptl_NItr sptl_Rglrztn;
global sptl_MskBrn sptl_MskObj;

global BCH fmriTEST fmriDIR; 

bboxes  = [   -78 78 -112 76  -50 85         
	      -64 64 -104 68  -28 72         
	      -90 91 -126 91  -72 109
              -95 95 -112 76 -50 95];
bbprompt =  [' -78:78 -112:76  -50:85  (Default)|'...
	     ' -64:64 -104:68  -28:72  (SPM95)   |'...
	     ' -90:91 -126:91  -72:109 (Template)|'...
             ' -95:95 -112:76  -50:95 '];
voxdims    = [ 1   1   1 ; 1.5 1.5 1.5 ; 2   2   2 ; 3   3   3 ; 4   4   4 ; 1   1   2 ; 2   2   4];
voxprompts = ' 1   1   1 | 1.5 1.5 1.5 | 2   2   2 | 3   3   3 | 4   4   4 | 1   1   2 | 2   2   4';
bases =     [0 0 0;2 2 2;2 3 2;3 3 3;3 4 3;4 4 4;4 5 4;5 5 5;5 6 5;6 6 6;6 7 6;7 7 7;7 8 7;8 8 8];
basprompt = 'none|2x2x2|2x3x2|3x3x3|3x4x3|4x4x4|4x5x4|5x5x5|5x6x5|6x6x6|6x7x6|7x7x7|7x8x7|8x8x8';


% fmri init stuff
%----------------
if nargin == 17
   rmflag = '';
elseif nargin == 18
   %ok
else
   error('This script needs 17 arguments (+ 1 optional rmflag)');   
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
if f99_exist([fmriDIR filesep 'LOG'],'NORMALIZE.OK')
   str = ['!rm -f ' fmriDIR filesep 'LOG' filesep 'NORMALIZE.OK'];
   eval(str);
end

% Start logging
%--------------
logfile = [fmriDIR filesep 'LOG' filesep 'normalize.log'];
tijd=spm('time');
disp('');
disp('*****************************************************************');
f99_log(logfile,['Spatial Normalization Script started at ' mat2str(tijd)]);
str = ['Start logging in ' logfile];
f99_log(logfile,str);
%--------------


% Setting GLOBAL values that are normaly set with 'Edit Defaults'
% ---------------------------------------------------------------
% customisation to yes:
sptl_CO = 1;
% starting estimates
sptl_Ornt = ase;
% ---------------------------------------------------------------


%if (nargin == 0)
% With no arguments, act as spm_sn3d_ui
SPMid = spm('FnBanner',mfilename,'2.20');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Normalize');
spm_help('!ContextHelp','spm_sn3d.m');

%-----------------------------------------------------------------------
%a1 = spm_input('Which option?',1,'m','Determine Parameters ...
%    Only|Write Normalised Only|Determine Parameters & Write Normalised',[1 2 3],3,'batch',{},'option');
a1 = option;

%nsubjects = spm_input('# Subjects','+1','w',1, ...
%    1,'batch',{},'nbsubjects');
nsubjects =1;
if nsubjects < 1,
    spm_figure('Clear','Interactive');
    return;
end;

object_masking = sptl_MskObj;
if (a1 == 1 | a1 == 3) %& sptl_CO ~= 1,
%       object_masking = spm_input('Mask object brain when registering?', '+1', 'm',...
%           'Mask object|Dont Mask object',[1 0],find([1 0] == object_masking),...
%           'batch',{},'object_masking');	
        object_masking = mask_object;
end;

% Select images..
%-----------------------------------------------------------------------
%for i=1:nsubjects,
i=1; % only implemented for 1 subject

subj(i) = struct('P','','PP','','objmask','','matname','');
  if a1 == 1 | a1 == 3,
%      if isempty(BCH),
%         subj(i).P = spm_get(1,'.img',['subj ' num2str(i) ...
%            ' - Image to determine parameters from']);
%      else,
%         subj(i).P = spm_input('batch', {},'image',i);
%      end;
      subj(i).P = fullfile(fmriDIR,objima);

      % object weight
      if object_masking==1,
%         if isempty(BCH),
%              subj(i).objmask = spm_get([0,1],'.img',...
%                ['Select object masking image (or Done for none)']);
%         else,
%              subj(i).objmask = spm_input('batch', {},objmask',i);
%         end;
         if isempty(maskima)
             subj(i).objmask = '';
         else
             subj(i).objmask = fullfile(fmriDIR,maskima);
         end;
      end;
      
      subj(i).matname = [spm_str_manip(subj(i).P,'sd') '_sn3d.mat'];
  else,
%      if isempty(BCH),
%         subj(i).matname = spm_get(1,'_sn3d.mat',['subj ' num2str(i) ...
%              ' - Normalisation parameter set:']);
%      else,
%         subj(i).matname = spm_input('batch', {},'matname',i);
%      end;
       subj(i).P = fullfile(fmriDIR,objima);
       subj(i).matname = [spm_str_manip(subj(i).P,'sd') '_sn3d.mat'];
  end;
    
  if a1 == 2 | a1 == 3,
%      if isempty(BCH)     
%         subj(i).PP = spm_get(Inf,'.img',['subj ' num2str(i) ...
%              ' - Images to write normalised']);
%      else
%         subj(i).PP = spm_input('batch', {},'images',i);
%      end
      subj(i).PP = f99_P(DIR,EXP);
  end;
%end;

if a1 == 1 | a1 == 3,

  % Get template(s)
  ok = 0;
  while ~ok,
%      if isempty(BCH),
%            Template = spm_get(Inf,'.img',['Template image(s)'],...
%               fullfile(SWD,'templates'));
%      else,
%            Template = spm_input('batch', {},'template');
%      end;
      Template = templima;
      vv=spm_vol(Template);
      if prod(size(vv))==1,
            ok = 1;
      elseif prod(size(vv)) ~= 0,
            tmp1 = cat(1,vv.dim);
            tmp2 = reshape(cat(3,vv.mat),4*4,prod(size(vv)));
            if ~any(any(diff(tmp1(:,1:3)))) & ~any(any(diff(tmp2,1,2))),
                ok=1;
            end;
      end;
  end;
  
  nbasis     = sptl_NBss;
  iterations = sptl_NItr;
  rglrztn    = sptl_Rglrztn;
  
  sptl_MskBrn = mask_templ;
  brainmask = sptl_MskBrn;

  
  
  
    % if sptl_CO ~= 1,
    % Customise the normalisation
    %-----------------------------------------------------------------------
    %a2  = spm_input('Normalisation Type?', '+1', 'm',...
    %    'Default Normalisation|Custom Normalisation', [0 1],1, 'batch',{},'type');
    a2 = 1;
    
    if a2 == 1,

      % Nonlinear options - # basis functions
      %-----------------------------------------------------------------------
      if prod(size(nbasis)) == 3,
        tmp = find(bases(:,1) == nbasis(1) & bases(:,2) == nbasis(2) & ...
            bases(:,3) == nbasis(3));
        if isempty(tmp), tmp = size(bases,1)+1; end
      else,
        tmp = size(bases,1)+2;
      end;
      %nb = spm_input('# Nonlinear Basis Functions?','+1','m',...
      %    [basprompt '|Custom'],[1:size(bases,1) 0], tmp,...
      %    'batch',{},'nonlin_func_nb');
      nb = nbf;
      
      if (nb>0), nbasis = bases(nb,:);
      elseif nb == 0,
        if (prod(size(sptl_NBss)) ~= 3) sptl_NBss = [5 6 5]; end;
        tmp = sprintf('%d %d %d', sptl_NBss(1), sptl_NBss(2), sptl_NBss(3));
        NBss = spm_input('# Basis Functions (x y z)','+0', 'w',...
            tmp, 3,'batch',{},'func_nb');
        NBss = NBss(:)';
        nbasis = NBss(:)';
      else, nbasis = 0; end;
      
      if prod(nbasis==0),
        iterations = 0;
        rglrztn    = 0;
      else,
          % Nonlinear options - # iterations
          %-----------------------------------------------------------------------
          tmp2 = [1 3 5 8 12 16];
          tmp = find(iterations == tmp2);
          if isempty(tmp) tmp = length(tmp2); end
          %iterations = spm_input('# Nonlinear Iterations?','+1','m',...
          %    ['1  nonlinear iteration |3  nonlinear iterations|' ...
          %      '5  nonlinear iterations|8  nonlinear iterations|' ...
          %      '12 nonlinear iterations|16 nonlinear iterations'],...
          %    tmp2, tmp,{},'nonlin_ite_nb');
          iterations = nli;
              
          % Get the amount of regularization
          %-----------------------------------------------------------------------
          tmp2 = [1 0.1 0.01 0.001 0.0001];
          tmp = find(tmp2 == rglrztn);
          if isempty(tmp) tmp = length(tmp2); end;
          %rglrztn = spm_input('Nonlinear Regularization','+1','m',...
          %    ['Extremely Heavy regularization|Heavy regularization|'...
          %      'Medium regularization|Light regularization|'...
          %      'Very Light regularization'], tmp2, tmp,...
          %    'batch',{},'nonlin_regular');
          rglrztn = nlr;
      end;

      def_brainmask = fullfile(SWD,'apriori','brainmask.img');
      %tmp = ~isempty(sptl_MskBrn);
      %if tmp, tmp = tmp + 1 - strcmp(sptl_MskBrn,def_brainmask); end;
      %tmp = spm_input('Mask brain when registering?', '+1', 'm',...
      %    'No Brain Mask|Default Brain Mask|Specified Brain Mask',[0 1 2],...
      %    tmp+1, 'batch',{},'mask_brain');
      tmp = mask_templ;
      if ~tmp, brainmask = ''; end;
      if tmp==1, brainmask = def_brainmask; end;
      %if tmp==2, brainmask = spm_get(1,'*.img','Specify brain ...
      %  mask'); end;
      if tmp==2, brainmask = spec_brainmask; end;
    end;
end;
%end;


if a1 == 2 | a1 == 3,

   % Get interpolation method (for writing images)
   %-----------------------------------------------------------------------
   %Hold = spm_input('Interpolation Method?','+1','m',...
   %    ['Nearest Neighbour|Bilinear Interpolation|'...
   %      'Sinc Interpolation (9x9x9)'],[0 1 -9], 2,...
   %    'batch',{},'interp');
   Hold = intplt;

   % Get bounding box.
   %-----------------------------------------------------------------------
   %if prod(size(sptl_BB)) == 6, bb = sptl_BB;
   %else,
   %  ans = spm_input('Bounding Box?','+1','m',...
   %      [ bbprompt '|Customise'], [1:size(bboxes,1) 0], 1,...
   %      'batch',{},'bounding_box');
   %  if ans>0, bb=reshape(bboxes(ans,:),2,3);
   %  else,
   %    directions = 'XYZ';
   %    bb = zeros(2,1);
   %    for d=1:3,
   %      str = sprintf('%d %d', bboxes(1,d*2-1), bboxes(1,d*2));
   %      bb(:,d) = spm_input(['Bounding Box ' directions(d) ],...
   %          '+1', 'e', str,2,'batch',{},...
   %          sprintf('direction%d',d));
   %    end;
   %  end;
   %end;
   bb = reshape(bbt,2,3);
   
   
   % Get output voxel sizes.
   %-----------------------------------------------------------------------
   %if prod(size(sptl_Vx)) == 3, Vox = sptl_Vx;
   %else,
   %  ans = spm_input('Voxel Sizes?','+1','m',...
   %      [ voxprompts '|Customise'], [1:size(voxdims,1) 0],3,...
   %      'batch',{},'voxel_sizes');
     
   %  if ans>0, Vox = voxdims(ans,:);
   %  else,
   %    Vox = spm_input('Voxel Sizes ','+0', 'e', '2 2 2',3,...
   %        'batch',{},'voxel_sizes_custom');
   %    Vox = reshape(Vox,1,3);
   %  end;
   %end;
   Vox = voxt;
else,
   bb     = sptl_BB;
   Vox    = sptl_Vx;
end;

% Go and do the work
%-----------------------------------------------------------------------
spm('Pointer','Watch')
if a1 == 1 | a1 == 3,
   f99_log(logfile,'DETERMINING REGISTRATION parameters');
   f99_log(logfile,['  on object image ' mat2str(subj(i).P)]);
   if f99_exist(deblank(spm_str_manip(subj(i).P,'h')),deblank(spm_str_manip(subj(i).P,'t'))) == 0
      f99_log(logfile,['    WARNING ! image', subj(i).P, ' does not exist !']);
   end
   if mask_object == 1
      f99_log(logfile,['   object image will be masked by ',subj(i).objmask]);
      if f99_exist(deblank(spm_str_manip(subj(i).objmask,'h')),deblank(spm_str_manip(subj(i).objmask,'t'))) == 0
         f99_log(logfile,['    WARNING ! image', subj(i).objmask, ' does not exist !']);
      end
   else
      f99_log(logfile,['     not masking object image']);
   end
   f99_log(logfile,['  to template ' Template]);
   if sptl_MskBrn == 1
      f99_log(logfile,['     template masked with ', brainmask]);
   else
      f99_log(logfile,['     not masking template image']);
   end
   f99_log(logfile,['  parameters in ' mat2str(subj(i).matname)]);
   f99_log(logfile,['  affine type is ' affinet]);
   f99_log(logfile,['  starting estimates ' mat2str(ase)]);
   if nbasis == [0 0 0]
      f99_log(logfile,['  affine only registration']);
   else
      f99_log(logfile,['  nonlinear basis ' mat2str(nbasis)]);
      f99_log(logfile,['  nonlinear iterations ' mat2str(iterations)]);
      f99_log(logfile,['  nonlinear regularization value ' mat2str(rglrztn)]);
   end
    
   if fmriTEST == 0
     for i=1:length(subj),
       spm('FigName',['Normalize: working on subj ' num2str(i)],Finter,CmdLine);
       fprintf('\rSubject %d:  ', i);
     
       f99_do_sn3d(subj(i).P,subj(i).matname,bb,Vox,...
           [nbasis(:)' iterations 8 rglrztn],Template,...
           brainmask, subj(i).objmask,affinet);
     end;
   end
end;

if a1 == 2 | a1 == 3,
   f99_log(logfile,'WRITING NORMALISED');
   f99_log(logfile,['  parameters in ' subj(i).matname]);
   if f99_exist(deblank(spm_str_manip(subj(i).matname,'h')),deblank(spm_str_manip(subj(i).matname,'t'))) == 0
      f99_log(logfile,['    WARNING ! mat file does not exist (yet) !']);
   end
   f99_log(logfile,['  bounding box is ' mat2str(bb)]);
   f99_log(logfile,['  with voxel sizes ' mat2str(Vox)]);
   f99_log(logfile,['  interpolation method ' mat2str(Hold)]);
   f99_log(logfile,'  on image(s) ');
   for j = 1:size(EXP)
      P = f99_P(DIR(j,:),EXP(j,:));
      f99_log(logfile,['   ' mat2str(DIR(j,:)) filesep mat2str(EXP(j,:))]);
      f99_checkP(DIR(j,:),EXP(j,:),1);
   end
   
   if fmriTEST == 0
     for i=1:length(subj),
       spm('FigName',['Normalize: writing subj ' num2str(i)],Finter,CmdLine);
       fprintf('\rSubject %d: Writing Normalized..', i);
	try
       		spm_write_sn(subj(i).PP,subj(i).matname,bb,Vox,Hold);
 	catch
		EXP
		DIR
		subj(i)
		spm_write_sn(subj(i).PP,subj(i).matname,bb,Vox,Hold);
	end    
     end;
   end
end;

fprintf('\r%60s%s', ' ',sprintf('\b')*ones(1,60));
spm('FigName','Normalize: done',Finter,CmdLine);
spm('Pointer');

statusfile = [fmriDIR filesep 'LOG' filesep 'NORMALIZE.OK'];
f99_log(statusfile,'');
if ~isempty(rmflag)
   f99_rm(rmflag,DIR,EXP);
end

tijd=spm('time');
f99_log(logfile,['Spatial Normalisation Script ended at ' mat2str(tijd)]);
disp('-----------------------------------------------------------------');
disp('');



