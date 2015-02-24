function tor_glm_specmask(TOR,varargin)
  
% Explicitly specifies a mask of voxels for glm esimation.
%
% Tor Wager modified this 12/16/01 to take inputs from structure input arg.
%
% FORMAT glm_specmask    
% Asks for files and steps to be done.
%
% FORMAT glm_specmask(cfg,mask)
% cfg,mask - fullpaths to the corresponding SPMcfg.mat and bigmask
% Modifies and an SPMcfg_file with a bigmask_file and estimates 
% statistics from it, overriding an previously existing results files 
% in the SPMcfg_file directory.
%                         
% This program creates a gray matter mask with all voxels to be used 
% for model estimation and estimates the model, if requested. Various 
% messages and promts will appear in the console - please read and 
% select the desired options from there.
%
% The program uses several SPM99 routines to 1) create a bigmask 
% image containing all gray matter voxels,  2) modify an SPMcfg.mat 
% file so that it refers to the mask image created in step 
% one as an explicit definition of which voxels to be included in the 
% statistics, and 3) esimate the GLM model from the modified SPMcfg.mat.
% The three steps can be performed separately, or in the same run of
% this program.
% 
% Creating the bigmask image involves segmenting the gray and white
% matter from an inplane anatomical image, combining them to produce
% skull-stripped brain, then smoothing this brain with an isotropic 
% kernel of 20 mm, and then setting all voxels with values > 0.1 to 1, 
% or 0 otherwise. The resulting image can be used as a mask. If a 
% different procedure for creating a mask image is preferred, it can 
% be performed separately.
%
% --------------------------------------------------------------------
% @(#)glm_specmask.m    1.2  Kalina Christoff    2000-11-13


CWD = pwd;
global SWD

if nargin ~=1 & nargin ~=3, 
  error('please specify either 1 or 3 input args');
end  

if nargin ==1;

% 1. Create a mask image if requested
% =======================================
%in_create = input(['\n\nWould you like to create a bigmask image file? [y|n]' ...
%                                         ' : '],'s');

in_create = TOR.in_create;

if strcmpi(in_create,'y'),
    % segment out the gray from the inplane, if no *_seg.img already
         % ---------------------------------------------------------------
    %in_segment = input(['\n\nDo you already have seg1 and seg2 inplane images? ',... 
    %                                              '[y|n] : '],'s');
    in_segment = TOR.in_segment;   

      if strcmpi(in_segment,'n');
                  spm fmri
                  %PF = spm_get(1,'.img',['Select (normalized) inplane image to' ...
                  %                                ' segment']);

		  PF = TOR.PF;
		  PG = TOR.PG;
                  %PG = spm_get(1,'.img',['Select the modality of image (T1.img or T2.img)'],...
                  %                                       fullfile(SWD,'templates'));

                  [segWD,brain] = fileparts(PF);  
                  cd(segWD);
                  spm_segment(PF,PG);

                  SEG = char(spm_get('files',segWD,[brain,'_seg1.img']),...
                                     spm_get('files',segWD,[brain,'_seg2.img']));
                  
    else SEG = spm_get(2,'*_seg*.img','Select seg1 and seg2 images');
       [segWD,brain] = fileparts(deblank(SEG(1,:))); 
                 brain = brain(1:end-5);
                 cd(segWD);
         
         end

         % create the bigmask - 1) extract brain by combining gray and white
    % matter, then 2) smooth brain_ image and 3) truncate at > 0.1
         % --------------------------------
    % 1)
         mode = 1; % save Extracted Brain only
         warning off; 
         spm_xbrain(SEG,mode); 
         warning on;
         
         % 2) smooth brain_ & create a bigmask.img file
    P = spm_get('files',segWD,['brain_',brain,'.img']);
         Q = 'bigmask';
         S = [20 20 20];
         spm_smooth(P,Q,S)
         
         
         % 3) trucate
         spm_imcalc_ui('bigmask.img','bigmask.img','i1>0.1');
         
         fprintf(['\n\nDone. An image "bigmask.img" was created in %s.\n', ...
                                 'It can be used as a mask at the statistics stage for this ',...
                                 'subject.\n\n'],segWD);
    bigmaskimg = spm_get('files',segWD,'bigmask.img');
         cd(CWD);

end

end


% 2. Modify the SPMcfg.mat file, if requested
% =============================================

if nargin == 3, 
  cfg_file  =varargin{1}; 
  bigmaskimg=varargin{2}; 
  in_modify='y';

elseif nargin ==1
  fprintf('\nIf you already have SPMcfg.mat file for this subject, you\n');
  fprintf('can modify it now and estimate the statistics with the full mask.\n\n');
  in_modify = TOR.in_modify;	% ...
  % input('Would you like to modify the SPMcfg.mat file? [y|n] : ','s');

end

if strcmpi(in_modify,'y'),
      if nargin == 1
                  %cfg_file = spm_get(1,'SPMcfg.mat',...
                  %                                                      'Select SPMcfg.mat file to modify',CWD);
		  cfg_file = TOR.cfg_file;

                  if ~exist('segWD','var');
                          %bigmaskimg=spm_get(1,'*.img',...
                          %                                         'Select a mask image with voxels to include');
				  bigmaskimg = TOR.bigmaskimg;
                  end
                end

           cfgWD=fileparts(cfg_file);
                segWD=fileparts(bigmaskimg);

                cd(cfgWD); load(cfg_file);
                xM.TH = ones(size(xM.TH))*(-Inf);
                vol=spm_vol([segWD,'/bigmask.img']);
           xM.VM=vol;
                
                save SPMcfg xM -append;
                fprintf('\n\nDone. The SPMcfg.mat file in \n%s\n',cfgWD);
                fprintf('has been modified. Estimating the model from it will\n');
                fprintf('include all voxels in the brain, as defined by the\n');
                fprintf('mask image file.\n\n');
end

                
% 3. Do the model estimation, if requested
% ========================================
if nargin == 1
  fprintf('\n\nYou can select to estimate the model now.\n');
  fprintf('Any previous results files in the SPMcfg.mat file directory\n');
  fprintf('will be overwritten.\n\n');
  %in_estim = ...
  %       input('Would you like to estimate the model now? [y|n] : ','s');

   in_estim = TOR.in_estim;

  if strcmpi(in_estim,'y') 
         if ~exist('cfg_file','var');
                cfg_file=spm_get(1,'SPMcfg.mat','Select SPMcfg.mat file to estimate');
      cfgWD = fileparts(cfg_file);
                cd(cfgWD); load(cfg_file);      
         end
    spm_spm(VY,xX,xM,F_iX0,Sess,xsDes)
  end
  
elseif nargin == 2
  
  spm_spm(VY,xX,xM,F_iX0,Sess,xsDes)

end
         








