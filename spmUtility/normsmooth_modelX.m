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
fmriDIR  = '/data/biman5'; 



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

if ~(exist('canonicalT1') == 1) | ~(exist('canonicalfunct') == 1) | ~(exist('segtemplate') == 1)
	% for jonides7
	% possible normalizeation template(s) - choose one of these.
	% -------------------------------------------------------------------
	canonicalT1 = '/data3/biman/spm99/spm99/templates/T1.img';
	canonicalfunct = '/data3/biman/spm99/spm99/templates/EPI.img';
	segtemplate = '/data3/biman/spm99/spm99/templates/T1_seg1.img';

	% for picasso
	% possible normalizeation template(s) - choose one of these.
	% -------------------------------------------------------------------
	canonicalT1 = '/usr/private/spm99/templates/T1.img';
	canonicalfunct = '/usr/private/spm99/templates/EPI.img';
	segtemplate = '/usr/private/spm99/templates/T1_seg1.img';
end

% bounding boxes - pick one in function call later.
defaultbb = [-78 78 -112 76  -50 85];
myotherbb = [-78 78 -150 56  -50 85];

% starting orientations - pick one in function call later.
% if you use -R in reconstruction, you're in Neuro - as of 8/23/01
neuro = [0 0 0 0 0 0 1 1 1 0 0 0];
radio = [0 0 0 0 0 0 -1 1 1 0 0 0];

% --------------------------------------------------------------------
% * if these don't exist, prompt...
% --------------------------------------------------------------------
%if ~(exist('t1name') == 1),
%	disp('Choosing the base image to calculate norm. parameters from.') 
%	t1name = input('Enter original t1 name, no homo corr, with .img extension','s');
%end

if ~(exist('substyle') == 1), 
	substyle = input('Enter subj index style, 1 for numbers, 2 for alphanumeric codes');
end

if ~(exist('subs') == 1), 
	error('You must enter subject numbers in subs variable before running.')
end

% --------------------------------------------------------------------
% * USER OPTIONS - enter here
% --------------------------------------------------------------------

Opt.wnorm = 'normal';
							% choices: 'segmented', 'normal'	for seg or normal
Opt.compute = 0;
Opt.allimgs = 1;
Opt.dort1 = 0;
Opt.fmap = 0;
Opt.tmap = 0;
Opt.con = 0;
Opt.bb = defaultbb;
Opt.orient = neuro;
Opt.wnorm = 'normal';
Opt.voxsize = [3.125 3.125 7];		% voxel size to write out
Opt.kernel = [12 12 12];			% smoothing kernel fwhm - should be 3 x voxel size for Gaussian Random Fields.
Opt.subs = subs;
% Opt.t1name = t1name;
Opt.rt1 = 'anatomy/ht1spgr.img';	% image to compute parameters from, with path after subject id directory

Opt



% --------------------------------------------------------------------
% * assign variable values from user options
% --------------------------------------------------------------------
wnorm = Opt.wnorm;
voxsize = Opt.voxsize;
kernel = Opt.kernel;

% list of processes to run - 1 is do it, 0 is don't do it.
compute = Opt.compute;	% calculate normalization parameters
allimgs = Opt.allimgs;	% 1 if you're normalizing all images, 0 for statistical maps (apply to all images?)
dort1 = Opt.dort1;		% apply to rt1? (or whichever anatomical you specify)
fmap = Opt.fmap;	% apply to f-maps?  these can be 0 if you only want nravols.
tmap = Opt.tmap;	% apply to t-maps?
con = Opt.con;	% apply to contrast maps?

% 1 for all images, 0 for tmaps/etc only.

% 03/08/01 Tor - I've been having some trouble using option 3, determine and
% 	write params.  It seems to work better if I use the rt1 image to determine
%	parameters from, and then apply them to tmaps, contast, and the rt1 itself
%  in a separate step.

% notes: remember to put in correct number of scans if you're normalizing all images!
%	if you don't, it will show up as an error getting file info on '' in spm_vol.



for JJ = subs
   if substyle == 2
   		% if subs get letter/number code   
   		snum = JJ{1};
	elseif substyle == 1
   		% if subs are numbered
   		snum = ['sub' num2str(JJ)];
   	else
		error('Unknown substyle - choose 1 or 2')
	end
   
   disp(['* Subject ' snum ' ______________________________________________________________'])

   % --------------------------------------------------------------------
   % * more normalization options - individual subject options
   % --------------------------------------------------------------------

   % possible images to determine params from
	% rt1 = [snum '/anatomy/rh' t1name];   % homocor,coreg with reslice, so this is realigned
	% rt1 = [snum '/anatomy/h' t1name];	   % raw homocor t1 name, no coregistration. does it work better?
	% segt1 = [snum '/anatomy/rh' t1name(1:end-4) '_seg1.img'];


	rt1 = fullfile(snum,Opt.rt1);
  
	% name of original, unsegmented t1 for application of normalization
	[rt1dir,rt1name,rt1ext] = fileparts(rt1);
	rt1name = [rt1name rt1ext];

   
   
   % directory where t, F, and con* images are
   if fmap ==1 | tmap == 1 | con == 1
	if ~exist('model') == 1, model = input('Enter analysis model number where con/t/F maps are:');,end
   	subdir = ['RESULTS/model' model filesep snum];
   end
   
   % actually use these images in the functions
   % ------------------------------------------
   if strcmp(wnorm,'segmented'), mytemplate = segtemplate;, elseif strcmp(wnorm,'normal'), mytemplate = canonicalT1;, 
   else error('incorrect wnorm'), end

   % img:										% compute parameters using this image
   if strcmp(wnorm,'segmented'), img = segt1;, elseif strcmp(wnorm,'normal'), img = rt1;, 
   else error('incorrect wnorm'), end

   dir = rt1dir;									% where to find the anatomicals
   exp = rt1name;									% the subject image you want to apply to (usually unsegmented).
   meanadir = [snum filesep 'task'];				% the directory where the scan1 subd's are, rel. to fmriDIR
 
   disp(['Template: ' mytemplate])
   disp(['Img to normalize: ' img])
   disp(['Apply norm to: ' exp])

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
      0,'',... 																					% m	asking
      [meanadir '/scan1';meanadir '/scan2';meanadir '/scan3';meanadir '/scan4';meanadir '/scan5';meanadir '/scan6';meanadir '/scan7';meanadir '/scan8';meanadir '/scan9'] ,... 																		% directories
         ['ravol*img';'ravol*img';'ravol*img';'ravol*img';'ravol*img';'ravol*img';'ravol*img';'ravol*img';'ravol*img'],... 																	% EXP: file name wildcards
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
            [meanadir '/scan1';meanadir '/scan2';meanadir '/scan3';meanadir '/scan4';meanadir '/scan5';meanadir '/scan6';meanadir '/scan7';meanadir '/scan8';meanadir '/scan9'] , ...
      ['nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img';'nravol*.img']) 
   
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
