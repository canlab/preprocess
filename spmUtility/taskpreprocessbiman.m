
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Obligatory first part of the automation script                           % 
%*************************************************************************** 
%  fmriTEST = 1 to run the simulation mode 
%           = 0 to run the actual post-processing 
%  fmriDIR  = the directory in which all your data is situated 
%--------------------------------------------------------------------------- 
% Last Modified: 03/14/01 by Tor Wager
%clear
global fmriTEST fmriDIR;    %set global variabels 
fmriTEST = 0; 

%subj = 12;  % enter this in workspace before running
ss = [7];
datadrive = [];

for subj = ss


if (exist('subj', 'var') ~= 1),subj = input('Enter subject number: ');,end
if (exist('datadrive', 'var') ~= 1),datadrive = input('Enter data drive for this subject: ');,end

subj = num2str(subj);
datadrive = num2str(datadrive);

fmriDIR  = (['/data' datadrive '/biman4/sub' subj]);

str = ['!mkdir ' fmriDIR '/task']
eval(str) 
str = ['cd ' fmriDIR]
eval(str)

!mv * task/
cd task

!chmod 777 *   % added because things were created with permissions spm couldn't read.


% Parameters
% ==================================================================
% Do these:
  recon = 1;
  slicetime = 1;
  realign = 1;
  smoothing = 0;
  moveem = 1;
 
TR = 2500;
TA = 2500;

% As of 4/16/01 fMRI_T = 35, fMRI_T0 = 17.
% Remember to set fMRI_T to # of slices in defaults
% and set fMRI_T0 to ref slice #.
% OR just set fMRI_T0 to fMRI_T * (refslice/numslices)
refslice = 17; 				% changed this
realign_mask = 0;			% changed this, too.
numscans = 6;
volsperscan = 272;
smooth_dims = [6 6 8];
moveSet = 'ravol';		% example: 'ravol', 'avol' to move this set of files.

% ==================================================================


% Commands
% ==================================================================

% Recon
if recon

!ls *.7 > file.txt
name = textread('file.txt','%s');
%name = name{1};  % for single scan recon only

% if no pfiles, attempt to unpack them.
%  if strcmp(name1(end-12:end),'pfiles.tar.gz')
if isempty(name)
disp(['Untarring/zipping pfiles.'])
!tar zxvf *.tar.gz
!ls *.7 > file.txt
name = textread('file.txt','%s');
name = name{4};

else
name = name{4};

end

str = ['!/data/gsp18.linux -m -v -R -n 64 ' name];
%str = ['!/data/gsp18b.linux -m -v -fx -fy -l -n 64 ' name];
disp(str); disp('Running....')
eval(str)
str = ['!/data/gsp18.linux -h -A -v -R -n 64 *.7'];		% outputs images in neuro format - R is R
%str = ['!/data/gsp18b.linux -h -A -v -fx -fy -l -n 64 *.7'];	
disp(str); disp('Running....')
eval(str)

% remove p-files afterward.  leave tarzipped pfiles.
!rm *.7

end

% Preprocessing
f99_CWD('LOG');             %make the LOG and RESULTS directory 
if slicetime
	f99_slice_timing(['task'],['vol*.img'],1,[],refslice,TR,TA);
        !gzip vol*img
end
if realign
	f99_realign(3,['task'],['avol*.img'],-9,3,realign_mask,0,1.00);
        !gzip avol*img
end
if smoothing
	f99_smooth(smooth_dims,'task','ravol*.img');
end

if moveem

% Copying images to scan directories
% copyimgs2scanTESTdirs(fmriDIR);				% sub-function

   str = ['cd ' fmriDIR '/task'];
   eval(str)
   
   [PI,nimgs] = getfiles([moveSet '*img']);
   % [PH,nimgs] = getfiles('*vol*hdr'); old slower code.
   % [PM,nimgs] = getfiles('*vol*mat');
   BN = PI{1}(1:end-8);
  disp(['Moving images for subject ' num2str(subj) ' to scan* dirs: basename is ' BN])

   index = 1;

   for KK = 1:numscans
      disp(['	Moving scan ' num2str(KK) '- images ' num2str(index) ' to ' num2str(index+volsperscan-1)])
      str = ['!mkdir scan' num2str(KK)'];
      eval(str)
      for JJ = index:index+volsperscan-1
	  
	  % set pad for filename
	  if JJ < 10, myzeros = '000';
	  elseif JJ < 100, myzeros = '00';
	  elseif JJ < 1000, myzeros = '0';
	  else myzeros = '';
	  end

          name = [BN myzeros num2str(JJ) '.*'];

          eval(['!mv ' name ' scan' num2str(KK) '/'])
          
	  %name = PH{JJ};
          %eval(['!mv ' name ' scan' num2str(KK) '/'])
          %name = PM{JJ};
          %eval(['!mv ' name ' scan' num2str(KK) '/'])
      end
      index = index+volsperscan;
   end
   

end

disp('Done!!!')

end % end loop through ss



%f99_coregister(1,0,... 
%    2,'Anatomy/t1.img',... 
%    5,'task/*mean*img',... 
%    ['task'],['avol*img']) 


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


% 3,... 
%            ['run01';'run02';'run03';'run04';... 
%             'run05';'run06';'run07';'run08'],... 
%            ['aim*.img';'aim*.img';'aim*.img';... 
%             'aim*.img';'aim*.img';'aim*.img';... 
%             'aim*.img';'aim*.img'],... 
%            -9,3,1,0,1.00,... 
%            '') 

% f99_coregister(1,0,...  
%    2,'Anatomy/t1.img',...
%    5,'task/meanavol_e3389_11_03_100_0001.img',...
%    ['task'],['ravol*img'])
  
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


% [5],... 
%           ['run01';'run02';'run03';'run04';... 
%            'run05';'run06';'run07';'run08'],... 
%           ['nra*.img';'nra*.img';'nra*.img';'nra*.img';... 
%            'nra*.img';'nra*.img';'nra*.img';'nra*.img'],... 
%           'SMOOTH.OK') 
  
  










