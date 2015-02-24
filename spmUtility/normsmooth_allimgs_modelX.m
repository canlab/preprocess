%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Obligatory first part of the automation script                           % 
%*************************************************************************** 
%  fmriTEST = 1 to run the simulation mode 
%           = 0 to run the actual post-processing 
%  fmriDIR  = the directory in which all your data is situated 
%--------------------------------------------------------------------------- 
global fmriTEST fmriDIR;    %set global variabels 
fmriTEST = 0; 
fmriDIR  = '/data4/intext2'; 
f99_CWD('LOG');             %make the LOG and RESULTS directory 


% 03/08/01 Tor - I've been having some trouble using option 3, determine and
% 	write params.  It seems to work better if I use the rt1 image to determine
%	parameters from, and then apply them to tmaps, contast, and the rt1 itself
%  in a separate step.

ss = [8 9 10 11 12];
drive = [2 2 4 4 4];

for JJ = ss
   snum = num2str(JJ);
   mydrive = num2str(drive(ss == JJ));
   
   try
      
      P = getfiles([fmriDIR filesep 'sub' snum filesep 'task' filesep 'nra*img']);
      if isempty(P)		% if no normalized files in there, start normalizing
         
         
         
         	% anatomical - determine parameters only
   			f99_sn3d(3,...																					% option
      		['../../data' mydrive '/intext/sub' snum '/Anatomy/rt1.img'], ...
      		0,'',... 																					% m	asking
      		['../../data' mydrive '/intext/sub' snum '/Anatomy'],... 						% directories
         ['rt1.img'],... 
         '/data3/biman/spm99/spm99/templates/T1.img', ...								% determine params from
         0,[],... 																				% masking of canonical				
         'affine3',[0 0 0 0 0 0 1 1 1 0 0 0],13,16,0.01,... 							% type and bounding box - neurological
         1,[-78 78 -122 76 -50 85],[4 4 4],'') 											% other params & voxel size - bilinear interp

  			% NOW APPLIES PARAMS FROM CANONICAL BRAIN AT 4 4 4 mm. smoothing 8 8 8. 
     
    		% apply to ravols 
    		% ================================================================================================================
    		f99_sn3d(2,...																						% option
      		['../../data' mydrive '/intext/sub' snum '/Anatomy/rt1.img'], ...
      		0,'',... 																							% m	asking
      		['sub' snum '/task'], ['ravol*img'], ...
      		'/data3/biman/spm99/spm99/templates/T1.img', ...											% determine params from
         0,[],... 																						% masking of canonical				
         'affine3',[0 0 0 0 0 0 1 1 1 0 0 0],13,16,0.01,... 									% type and bounding box
         1,[-78 78 -122 76 -50 85],[4 4 4],'') 													% other params & voxel size
      
      
      		%['sub' snum '/task/scan1'; ...
       	%  'sub' snum '/task/scan2';...
        	%'sub' snum '/task/scan4';...        
			%'sub' snum '/task/scan5';...       
			%'sub' snum '/task/scan6';...      
			%'sub' snum '/task/scan7';...       
			%'sub' snum '/task/scan8'...
      		%],...  % directories

        	% ['ravol*img';'ravol*img';'ravol*img';'ravol*img';'ravol*img';'ravol*img';'ravol*img';'ravol*img'],... 
      
      		% apply to anatomical
    		% ================================================================================================================
    		f99_sn3d(3,...																					% option
      		['../../data' mydrive '/intext/sub' snum '/Anatomy/rt1.img'], ...
      		0,'',... 																					% m	asking
      		['../../data' mydrive '/intext/sub' snum '/Anatomy'],...						% directories
         ['rt1.img'],... 
         '/data3/biman/spm99/spm99/templates/T1.img', ...								% determine params from
         0,[],... 																				% masking of canonical				
         'affine3',[0 0 0 0 0 0 1 1 1 0 0 0],13,16,0.01,... 							% type and bounding box - neurological
         1,[-78 78 -122 76 -50 85],[4 4 4],'') 											% other params & voxel size - bilinear interp
      

      
         Q = getfiles([fmriDIR filesep 'sub' snum filesep 'task' filesep 'ravol*gz']);
      		if isempty(A)
				disp('Zipping RAVOLS')
				eval(['!gzip ' fmriDIR filesep 'sub' snum filesep 'task' filesep 'ravol*'])            
      		else 
         		disp(['Zipped ravol files already exist!  Skipping compression...'])
         end
            
   else
   		disp(['Normalized files already exist!  Skipping normalization...'])
   end % end if isempty(P)
   
   
   
   
   % smoothing   
   % ================================================================================================================  
   try
      P = getfiles([fmriDIR filesep 'sub' snum filesep 'task' filesep 'snra*img']);
      if isempty(P)		% if no snra files in there, start smoothing

         f99_smooth([8 8 8],	['sub' snum '/task'], ['nravol*img']);
         
      else 
         disp(['Smoothed nravol files already exist!  Skipping smoothing...'])
         
      end
      
      Q = getfiles([fmriDIR filesep 'sub' snum filesep 'task' filesep 'snra*gz']);
      if isempty(Q)	
         disp('Zipping NRAVOLS')
         eval(['!gzip ' fmriDIR filesep 'sub' snum filesep 'task' filesep 'nravol*']) 	
         
      else
         disp(['Zipped nravol files already exist!  Skipping compression...'])
         
      end
      
	catch
      warning('Problem Smoothing!')
   end
   
catch
   disp('Problem normalizing!')
end   


end % loop thru subjects



   
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
