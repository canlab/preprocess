%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Obligatory first part of the automation script                           % 
%*************************************************************************** 
%  fmriTEST = 1 to run the simulation mode 
%           = 0 to run the actual post-processing 
%  fmriDIR  = the directory in which all your data is situated 
%--------------------------------------------------------------------------- 
global fmriTEST fmriDIR;    %set global variabels 
fmriTEST = 1; 
fmriDIR  = '/images_rinda3/beatse/sidev6'; 
f99_CWD('LOG');             %make the LOG and RESULTS directory 


%*************************************************************************** 
% 'f99_slice_timing(DIR,EXP,Seq,sl_order,refslice,timing)' 
%*************************************************************************** 
%  'DIR'       - [matrix containing directories of experiments] 
%  'EXP'       - [matrix containing directories of images wildcard] 
%  'Seq'       - slice acquisition order (1,2,3 = asc, desc, interl,user 
%                                         specified) 
%  'sl_order'  - chronological order of slice acquisition (if Seq = 4) 
%  'refslice'  - slice for time 0 
%  'TR'        - repitition time 
%  'TA'        - acquisition time (time between start of acq of first 
%                slice and the end of acq of the last slice) 
%====================================================================== 
f99_slice_timing(... 
     ['run01';'run02';'run03';'run04';... 
      'run05';'run06';'run07';'run08'],... 
     ['im*.img';'im*.img';... 
      'im*.img';'im*.img';... 
      'im*.img';'im*.img';... 
      'im*.img';'im*.img'],... 
     1,[],... 
     1,... 
     3500,2880) 
  

%*************************************************************************** 
% 'f99_coregister(option, target_mod, targetima, object_mod,...'            % 
%                'objectima, DIR, EXP)'                                     % 
%*************************************************************************** 
% COREGISTER 2 images 
% 'do what?' 
%    1 - coregister only 
%    2 - reslice only 
%    3 - coregister & reslice 
% 'Use Mutual Information?' 
%    1 = yes ; 0 = no 
% 'modality of target image' 
%    1 = PET    | 2 = T1 MRI | 3 = T2 MRI 
%    4 = PD MRI | 5=  EPI    | 6 = Transm 
% 'target image' 
%    e.g. 'run1/first_epi.img' 
% 'modality of object image' 
%    1 = PET    | 2 = T1 MRI | 3 = T2 MRI 
%    4 = PD MRI | 5=  EPI    | 6 = Transm 
% 'object image' 
%    e.g. 'anat/mpr.img' 
% 'Other images to reslice' 
% DIR   - [matrix containing directories of experiments0] 
% EXP   - [matrix containing directories of images wildcard] 
%============================================================================= 
f99_coregister(1,0,... 
    5,'run01/aim_001.img',... 
    2,'anato/mpr.img',... 
    [],[]) 
  

%*************************************************************************** 
% 'f99_realign(do_what,DIR,EXP,reslice_type,create_what,mask_images,...'   * 
%             'adj_se,real_qual,rmflag)'                                   * 
%*************************************************************************** 
%   do_what? 
%         1 = coregister only 
%         2 = reslice only 
%         3 = coregister and reslice 
%   DIR   - [matrix containing directories of experiments0] 
%   EXP   - [matrix containing directories of images wildcard] 
%   reslice_type? 
%         1  = Trilinear 
%        -9  = Sinc 
%        Inf = Fourrier (only if isotropic voxels) 
%   create_what? 
%         1 = All Images (1..n)       | 2 = Images 2..n 
%         3 = All Images + Mean Image | 4 = Mean Image Only' 
%   mask_images?            -  1 = yes, 0 = no 
%   adj_sampling_errors?    -  1 = yes, 0 = no 
%   registration quality?   -  1.00 to 0.001 
%         'Quality 1.00  (slowest/most accurate)  |Quality 0.90|' ... 
%         'Quality 0.75 |Quality 0.50|Quality 0.25|Quality 0.10|' ... 
%         'Quality 0.05 |Quality 0.01|' ... 
%         'Quality 0.005|Quality 0.001 (fastest/poorest)' 
%---OPTIONALY-------------------------------------------------------------- 
%   rmflags - you can specify 'REALIGN.OK' as last parameter to delete the 
%             original images after completion of realigment - USE WITH CAUTION 
%============================================================================= 
f99_realign(3,... 
            ['run01';'run02';'run03';'run04';... 
             'run05';'run06';'run07';'run08'],... 
            ['aim*.img';'aim*.img';'aim*.img';... 
             'aim*.img';'aim*.img';'aim*.img';... 
             'aim*.img';'aim*.img'],... 
            -9,3,1,0,1.00,... 
            '') 

  
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
f99_sn3d(3,...
         'anato/mpr.img',0,'',... 
         ['run01';'run02';'run03';'run04';... 
          'run05';'run06';'run07';'run08'],... 
         ['ra*.img';'ra*.img';'ra*.img';'ra*.img';... 
          'ra*.img';'ra*.img';'ra*.img';'ra*.img'],... 
         '/usr/local_rinda/src/SPM/spm99b/templates/T1.img',0,[],... 
         'affine3',[0 0 0 0 0 0 -1 1 1 0 0 0],13,16,0.01,... 
         1,[-78 78 -122 76 -50 85],[2 2 2],'') 
  

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
f99_smooth([5],... 
           ['run01';'run02';'run03';'run04';... 
            'run05';'run06';'run07';'run08'],... 
           ['nra*.img';'nra*.img';'nra*.img';'nra*.img';... 
            'nra*.img';'nra*.img';'nra*.img';'nra*.img'],... 
           'SMOOTH.OK') 
  
  