%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Obligatory first part of the automation script                           % 
%*************************************************************************** 
%  fmriTEST = 1 to run the simulation mode 
%           = 0 to run the actual post-processing 
%  fmriDIR  = the directory in which all your data is situated 
%--------------------------------------------------------------------------- 
global fmriTEST fmriDIR;    %set global variabels 
fmriTEST = 0; 

% for jonides7
fmriDIR  = '/data4/intext2'; 

% for picasso
fmriDIR  = '/data/biman4'; 
fmriDIR = '/data2/inhib';


f99_CWD('LOG');             %make the LOG and RESULTS directory 

eval(['cd ' fmriDIR])

% --------------------------------------------------------------------
% * some options for images to use, etc. are listed here.
% * check below to see what is actually used.
% --------------------------------------------------------------------

% notes:
%	7/2/01: rt1 -> T1.img produces very bad results
%				meanavol -> T1.img looks ok with check reg.
%				meanavol -> EPI.img looks ok except for ss 11 and 12 - lateral moves.

% for jonides7
% possible normalizeation template(s) - choose one of these.
% -------------------------------------------------------------------
canonicalT1 = '/data3/biman/spm99/spm99/templates/T1.img';
myavg = '/data1/intext/avganat/avgt1.img';
canonicalfunct = '/data3/biman/spm99/spm99/templates/EPI.img';
segtemplate = '/data3/biman/spm99/spm99/templates/T1_seg1.img';

% for picasso
% possible normalizeation template(s) - choose one of these.
% -------------------------------------------------------------------
canonicalT1 = '/usr/private/spm99/templates/T1.img';
canonicalfunct = '/usr/private/spm99/templates/EPI.img';
segtemplate = '/usr/private/spm99/templates/T1_seg1.img';


% use this template
mytemplate = canonicalT1;

% bounding boxes - pick one in function call later.
defaultbb = [-78 78 -112 76  -50 85];
myotherbb = [-78 78 -150 56  -50 85];

% starting orientations - pick one in function call later.
% if you use -R in reconstruction, you're in Neuro - as of 8/23/01
neuro = [0 0 0 0 0 0 1 1 1 0 0 0];
radio = [0 0 0 0 0 0 -1 1 1 0 0 0];

% voxel size to write out
voxsize = [3.125 3.125 4];
voxsize = [2 2 2];

% smoothing kernel fwhm - should be 3 x voxel size for Gaussian Random Fields.
kernel = [12 12 12];
kernel = [8 8 8];

% list of processes to run - 1 is do it, 0 is don't do it.
compute = 0;	% calculate normalization parameters
allimgs = 1;	% 1 if you're normalizing all images, 0 for statistical maps (apply to all images?)
dort1 = 0;	% apply to rt1? (or whichever anatomical you specify)
fmap = 0;	% apply to f-maps?  these can be 0 if you only want nravols.
tmap = 0;	% apply to t-maps?
con = 0;	% apply to contrast maps?

% 1 for all images, 0 for tmaps/etc only.

% 03/08/01 Tor - I've been having some trouble using option 3, determine and
% 	write params.  It seems to work better if I use the rt1 image to determine
%	parameters from, and then apply them to tmaps, contast, and the rt1 itself
%  in a separate step.

% notes: remember to put in correct number of scans if you're normalizing all images!
%	if you don't, it will show up as an error getting file info on '' in spm_vol.

for JJ = subs
   snum = num2str(JJ);
   
   % --------------------------------------------------------------------
   % * more normalization options - individual subject options
   % --------------------------------------------------------------------

   % possible images to determine params from
	rt1 = ['sub' snum '/anatomy/rhft1.img'];   % homocor,coreg with reslice, so this is realigned
	segt1 = ['sub' snum '/anatomy/rhft1_seg1.img'];

   % name of original, unsegmented t1 for application of normalization
	[rt1dir,rt1name,rt1ext] = fileparts(rt1);
	rt1name = [rt1name rt1ext];

   
   
   myP = getfiles(['sub' snum '/task/meanavol*'])
   meana = myP{2};
   meanadir = ['sub' snum '/task'];
   % meanaexp = 'mean.img'; not used here - use when normalizing mean to template.
   
   % directory where t, F, and con* images are
   subdir = ['RESULTS/model' model '/sub' snum];
   
   % actually use these images in the functions
   % ------------------------------------------
	img = '/data2/inhib/011018cg/anatomy/homocor.img';
	dir = '/data2/inhib/011018cg/anatomy/';
	exp = [];
	targetdir = ['/data2/inhib/011018cg/flanker/scan5';'/data2/inhib/011018cg/flanker/scan6'];

   %img = segt1;										% compute parameters using this image
   %dir = rt1dir;									% where to find the anatomicals
   %exp = rt1name;									% the subject image you want to apply to.
   
   %if JJ > 7,fmriDIR = '/data2/intext';,f99_CWD('LOG'); ,end
	%if JJ > 9,fmriDIR = '/data4/intext';,f99_CWD('LOG'); ,end
   
   % --------------------------------------------------------------------
   
   
   % notes:
   % FIL recommends masking the template, but not the object
   % i've chosen 'mask template' with spm99 apriori brainmask (i assume).
   % might also use segmented gray matter from template and object...
   % but of course, you need to normalize in order to segment.
   % so this would be a kind of recursive process.
   %
   % smoothness of the con maps should be 3 x the voxel size, so
   % we're smoothing to 12 mm so that random field theory will work.
   % Tom doesn't smooth at all - he's resampled at 1/3 original vox size instead.
   % more smoothness OR more df = lower thresholds in gaussian RFT.
   
   if compute
   % anatomical - determine parameters only
   % ================================================================================================================
   f99_sn3d(1,...																				% option
      img, ...																					% determine params from
      0,'',... 																					% masking
      dir,... 																					% directories
         [exp],... 																				% EXP: file name wildcards
         mytemplate, ...																		% template
         1,[],... 																				% masking of canonical				
         'affine3',neuro,13,16,0.01,... 													% type and bounding box - neurological
         1,defaultbb,voxsize,'') 															% other params & voxel size - bilinear interp
   
 
   end
   
   
   if allimgs
         
      % apply to ravols
      % ================================================================================================================
    	f99_sn3d(2,...																				% option
      img, ...																					% determine params from
      0,'',... 																					% masking
      [targetdir(1,:);targetdir(2,:)] ,... 																				% directories
         ['ravol*img';'ravol*img'],... 																	% EXP: file name wildcards
         [], ...																				% template
         0,[],... 																				% masking of canonical				
         'affine3',neuro,13,16,0.01,... 													% type and bounding box - neurological
         1,defaultbb,voxsize,'') 															% other params & voxel size - bilinear interp
   
   end
      
   if dort1
      
% apply to rt1
% ================================================================================================================
   f99_sn3d(2,...																				% option
      img, ...																					% determine params from
      0,'',... 																					% m	asking
      rt1dir,... 																				% directories
         [exp],... 																				% EXP: file name wildcards
         [], ...																				% template
         0,[],... 																				% masking of canonical				
         'affine3',neuro,13,16,0.01,... 													% type and bounding box - neurological
         1,defaultbb,voxsize,'') 															% other params & voxel size - bilinear interp
      
   end
   

   if fmap
      
      % apply to F maps 
    % ================================================================================================================
    f99_sn3d(2,...																				% option
      img, ...																					% determine params from
      0,'',... 																					% m	asking
      subdir,... 																				% directories
         ['spmF*img'],... 																	% EXP: file name wildcards
         [], ...																				% template
         0,[],... 																				% masking of canonical				
         'affine3',neuro,13,16,0.01,... 													% type and bounding box - neurological
         1,defaultbb,voxsize,'') 															% other params & voxel size - bilinear interp
   
   end
   
   if tmap
      
      % apply to T maps 
    % ================================================================================================================
    f99_sn3d(2,...																				% option
      img, ...																					% determine params from
      0,'',... 																					% m	asking
      subdir,... 																				% directories
         ['spmT*img'],... 																	% EXP: file name wildcards
         [], ...																				% template
         0,[],... 																				% masking of canonical				
         'affine3',neuro,13,16,0.01,... 													% type and bounding box - neurological
         1,defaultbb,voxsize,'') 															% other params & voxel size - bilinear interp
   
   end
   
   if con
      
      % apply to contrast maps
      % ================================================================================================================
    f99_sn3d(2,...																				% option
      img, ...																					% determine params from
      0,'',... 																					% masking
      subdir,... 																				% directories
         ['con*img'],... 																		% EXP: file name wildcards
         [], ...																				% template
         0,[],... 																				% masking of canonical				
         'affine3',neuro,13,16,0.01,... 													% type and bounding box - neurological
         1,defaultbb,voxsize,'') 															% other params & voxel size - bilinear interp
          
   end
   
   if allimgs
      % smoothing   
   % ================================================================================================================  
   try
      f99_smooth(kernel,... 
            [meanadir '/scan1';meanadir '/scan2';meanadir '/scan3';meanadir '/scan4';meanadir '/scan5';meanadir '/scan6';meanadir '/scan7';meanadir '/scan8'] , ...
      ['nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img']) 
   
   		% str = ['!\rm ' meanadir '/scan*/nravol*img'] 
   		% eval(str)
   
	catch
   
	end

end

if fmap
   % smoothing   
   % ================================================================================================================  
   f99_smooth(kernel,... 
      ['RESULTS/model' model '/sub' snum], ...
      ['nspmF*.img']) 
end

if con
  f99_smooth(kernel,... 
      ['RESULTS/model' model '/sub' snum], ...
      ['ncon*.img']) 
end

if tmap
   f99_smooth(kernel,... 
      ['RESULTS/model' model '/sub' snum], ...
      ['nspmT*.img'])
end

     
end % loop through subjects



   
   %f99_sn3d(3,...																						% option
   %   ['sub' snum '/Anatomy/rt1.img'], ...
   %   0,'',... 																							% m	asking
   %   ['sub' snum '/Anatomy'],... 																	% directories
   %      ['rt1.img'],... 
   %      '/data3/biman/spm99/spm99/templates/T1.img', ...										% determine params from
   %      0,[],... 																						% masking of canonical				
   %      'affine3',[0 0 0 0 0 0 -1 1 1 0 0 0],13,16,0.01,... 								% type and bounding box
   %      -9,[-78 78 -122 76 -50 85],[3.125 3.125 4],'') 										% other params & voxel size
      
   % t-maps   
   % ================================================================================================================
   %f99_sn3d(3,...																						% option
   %   ['RESULTS/model' model '/sub' snum '/spmT_0003.img'], ...
   %   0,'',... 																							% m	asking
   %   ['RESULTS/model' model '/sub' snum],... 												% directories
   %      ['spmT*.img'],... 
   %      '/data3/biman/spm99/spm99/templates/T1.img', ...										% determine params from
   %      0,[],... 																						% masking of canonical				
   %      'affine3',[0 0 0 0 0 0 -1 1 1 0 0 0],13,16,0.01,... 								% type and bounding box
   %      -9,[-78 78 -122 76 -50 85],[3.125 3.125 4],'') 										% other params & voxel size

   
   
%*************************************************************************** 
% 'f99_sn3d(option,objima,mask_object,maskima,DIR,EXP,templima,...' 
%          'mask_templ,spec_brainmask,affinet,ase,nbf,nli,nlr,...' 
%          'intplt,bbt,voxt,rmflag)' 
%*************************************************************************** 
% 'Do What?' 
%        1 = 'Determine Parameters Only' 
%        2 = 'Write Normalised Only' 
%        3 = 'Determine Parameters & Write Normalised' 
% 'Image to determine parameters from (objectimage)', OR 
%        in case of write normalised only, image for which parameters exist 
%        e.g. 'anat/anatmy.img' 
% 'Mask the object?' 
%        1 = yes, 0 = no 
% 'Image that contains the mask for the objectimage' 
%        e.g. 'anat/noskull.img' 
% 'DIR'  - [matrix of dirs] 
% 'EXP'  - [matrix of experiments0] 
% 'template image'  -  [one or more templates] 
%        The image will be fitted (in the least 
%        squares sense) to the optimum linear combination of these templates. 
% 'Mask the template image?' 
%       0 = 'no masking' 
%       1 = 'mask the template image with the mask in SPM/apriori directory' 
%       2 = 'mask with the spec_brainmask (user defined)' 
% 'user specified brainmask for template' 
% 'Affine Normalisation Type?'  - see spm_affsub3.m for details on choices 
%        'rigid1' - rigid body normalisation 
%        'rigid2' 
%        'rigid3' 
%        'affine1' - affine normalisation 
%        'affine2' 
%        'affine3'  - this is default SPM99 
% 'Affine Starting Estimates?' 
%        'Neurological Convention (R is R)' 
%             The starting estimates do not include a left right flip. 
%             [0 0 0 0 0 0  1 1 1 0 0 0]. 
%        'Radiological Convention (L is R)' 
%             [0 0 0 0 0 0 -1 1 1 0 0 0]. 
% '# Basis Functions (x y z)' 
%        1=none  -  USE THIS IF AFFINE ONLY NORMALISATION WANTED 
%        2=2x2x2 | 3=2x3x2 | 4=3x3x3 | 5=3x4x3 | 6=4x4x4 | 7=4x5x4 | 8=5x5x5 
%        9=5x6x5 | 10=6x6x6 | 11=6x7x6| 12=7x7x7 | 13=7x8x7 | 14=8x8x8 
% '# Nonlinear Iterations?' 
%        can be set to : 1, 3, 5, 8, 12 or 16 
% 'Nonlinear Regularization' 
%        can be set to : 1, 0.1, 0.01, 0.001 or 0.0001 
%        from 'Extremely Heavy regularization' to 'Very Light regularization' 
%        0.01 is SPM99 default, use Heavy if you get distorted results0 
% 'Interpolation Method?' 
%         0 = Nearest Neighbour | 1 = Bilinear Interpolation |... 
%        -9 = Sinc Interpolation (9x9x9) 
% 'Bounding Box?' 
%        [-x x -y y -z z] in mm e.g. [-78 78 -112 76  -50 85]  (Default SPM99) 
% 'Voxel Sizes ' 
%        The voxel sizes (x, y & z, in mm) of the written normalised images. 
%        [x y z] in mm e.g. [2 2 2] 
%---OPTIONALLY--------------------------------------------------------------- 
%   rmflags - you can specify 'NORMALIZE.OK' as last parameter to delete the 
%                unsmoothed images after completion - USE WITH CAUTION 
%--------------------------------------------------------------------------- 

%******************************************************************************* 
%f99_smooth(SP,DIR,EXP,rmflag)                                                 % 
%******************************************************************************* 
% SP    - smoothing kernel parameters [x y z] in mm 
% DIR   - [matrix of directories containing the images] 
% EXP   - [matrix of wildcard used to select the images in the DIR directories] 
%---OPTIONALLY--------------------------------------------------------------- 
%   rmflags - you can specify 'SMOOTH.OK' as last parameter to delete the 
%                unsmoothed images after completion - USE WITH CAUTION 
%--------------------------------------------------------------------------- 
