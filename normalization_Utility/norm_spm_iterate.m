function norm_spm_iterate(imagenames,iter,template,subjnames,varargin)
    % norm_spm_iterate(imagenames,iter,template,subjnames)
    %
    % Iterative spm normalization to the mean image in a set
    % Start with UNnormalized images and normalization parameters from SPM
    % Then run this on them to improve warping.
    %
    % imagenames = filenames('*/T1.img','absolute','char');
    % prmnames = filenames('*/*_sn.mat','absolute','char');
    % norg = 10;
    % numgen = 10;
    % iter = 10;
    % norm_spm_iterate(imagenames,prmnames,norg,numgen,iter)
    %
    % Optional inputs
%                         case 'mask', mask = varargin{i+1}; varargin{i+1} = [];
% 
%                     case 'graphics', dographics = 1;
%                         
%                     case 'i_start', i_start=varargin{i+1};
    %
    % Mask is interpolated to the space defined by the first image
    % and saved as mask.img
    %
    % norm_spm_iterate(imagenames,prmnames,norg,numgen,iter,'mask',maskname)
    %
    % Tor Wager
    % Nov. 06, version 11/26
    %
    % Create a mask image from previous mean:
    % img = 'mean.img'
    % V = spm_vol(img);
    % v = spm_read_vols(V);
    % v = v > 10;
    % VO = V; VO.fname = 'startmask.img';
    % spm_write_vol(VO,v);

    
    diary Norm_spm_Log.txt

    % -----------------------------------------------------
    % SETUP mat files and mask image
    % -----------------------------------------------------
    % read param files, make sure we can load them,
    % save copies in norm_subj???_sn3d.mat
    % replace prmnames with these names


    %[dographics,mask,template,VG,maskInfo,n,matnames,VF] = setup_norm;
    setup_norm;

    % variables are passed out of inline fcns only if they're used in the
    % main function.
    % passing in dographics causes this var to be returned by setup_norm.
    setup_spm_norm_options(dographics);     % returns flags for spm normalization




    % -----------------------------------------------------
    % Iterate steps of mean image creation and warping
    % -----------------------------------------------------
    fitness_by_iteration = zeros(1,n);

    for iteration = i_start:iter

        fprintf(1,'_____________________________________________________\n\n')
        fprintf(1,'Starting Iteration : %3.0f\n', iteration)
        flags.smosrc = source_w(iteration);
        flags.cutoff = cutoff_w(iteration);
        fprintf(1,'Source smoothing: %3.2f, Spatial basis cutoff: %3.2f\n',flags.smosrc,flags.cutoff);
        fprintf(1,'_____________________________________________________\n\n')


        % -----------------------------------------------------
        % Apply normalization
        % -----------------------------------------------------
        %
        mkdir(['iteration' num2str(iteration)]);
        cd(['iteration' num2str(iteration)]);
        
        for s = 1:n
            
            % load .mat file name from last iteration and get params
            % if it > 1
            % -------------------------------------------------
            T = [];
%             if iteration > 1
                % load params for this subject
%                 myprmfile = deblank(matnames{s});
%                 load(myprmfile, 'T')
%             end
            
             % update .mat file name for this iteration
            % if it > 1
            % -------------------------------------------------
            [curdir,curfname,ext] = fileparts(matnames{s});
            matnames{s} = fullfile(pwd, [curfname ext]);
            
            % *** we have to make sure this writes out an .img and .mat file
            % if matnames is entered, will save sn3d.mat file in this file
            params = spm_normalise_mod(VG,VF(s),matnames{s},['..' filesep 'mask.img'],[],flags, T);

            
            % *** define warped (spgr) image names in CURRENT DIR.  
            % save both image and norm. mat file
            spm_write_sn(VF(s),params,wflags)
            [A,B]=fileparts(imagenames(s,:));
            movefile([A filesep 'w' B '.img'],[subjnames{s} '_wspgr.img']);
            movefile([A filesep 'w' B '.hdr'],[subjnames{s} '_wspgr.hdr']);
            
        end



        % Create mean image
        % -----------------------------------------------------
        % write mean.img, which is the new template
        %
        % also load param files for all subjects and return warping info for all
        % subjects in cells of prm.

        prmnames = str2mat(matnames{:});
        
        

        % this applies CURRENT params to ORIGINAL images and writes mask of trimmed mean
        [meandata,fitness,prm] = mean_warped_image(imagenames,prmnames,'mask',['..' filesep 'mask.img']);

        imagenames = filenames('*_wspgr.img','absolute','char');
        
        fitness_by_iteration(iteration,:) = fitness;

        
        % Update image and mat names
        % ----------------------------------------------------- 
        % *** (keeping same matnames will overwrite old mat files)
        % re-map vols for next iteration
        % VF and imagenames should always remain the same (the orig.
        % images)
        VF = spm_vol(imagenames);
        
        template = filenames('mean.img','absolute','char');
        
        %*****not sure if this is needed
%         fprintf(1,'Found template image: \n%s\n\n', template);
%         VG = spm_vol(template);
%         disp('Creating mask.img in same space as template');
%         fprintf(1,'Found input image: %s\n', mask);
%         maskdat = iimg_reslice(template,mask,'write','outname','mask.img');
        
        banner;
        cd ..

    end  % iterations


    % write final mean image and individual normalized images

    prmnames = str2mat(matnames{:});
    [meandata,fitness,prm] = mean_warped_image(imagenames,prmnames,'write');

    fitness_by_iteration(end+1,:) = fitness;

    diary off

    %%% end main function %%%






    %%% Inline functions %%%

    % ___________________________________________________________________
    %
    % Inline: Print banner for start of one iteration
    % ___________________________________________________________________

    function banner

        fprintf(1,'_____________________________________________________\n\n')
        fprintf(1,'Completed Iteration : %3.0f\n', iteration)

        fprintf(1,'Fitness by iteration: mean\n')
        fprintf(1,'%3.4f ',mean(fitness_by_iteration,2)');
        fprintf(1,'\n')

        fprintf(1,'Fitness by iteration: min\n')
        fprintf(1,'%3.4f ',min(fitness_by_iteration,[],2)');
        fprintf(1,'\n')

        fprintf(1,'\n_____________________________________________________\n')

    end




    % ___________________________________________________________________
    %
    % Inline: Setup mask, population type, template for Normalization
    % ___________________________________________________________________

    function setup_norm

        dographics = 0;
        i_start=1;

        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}
                    case 'mask', mask = varargin{i+1}; varargin{i+1} = [];

                    case 'graphics', dographics = 1;
                        
                    case 'i_start', i_start=varargin{i+1};

                    otherwise
                        error('Unknown string option');
                end
            end
        end

        % template image
        [dd,ff,ee] = fileparts(template);
        template = which([ff ee]);
        fprintf(1,'Found template image: \n%s\n\n', template);
        VG = spm_vol(template);
        %[volInfo,dat_template] = iimg_read_img(template,2);

        % mask image
        % make sure this is in the space defined by the first image
        disp('Creating mask.img in same space as template');
        if exist('mask','var')
            fprintf(1,'Found input image: %s\n', mask);
            maskdat = iimg_reslice(template,mask,'write','outname','mask.img');
        else
            fprintf(1,'Mask is all ones\n');
            maskdat = ones(VG.dim(1:3));
            VO = VG;
            VO.fname = 'mask.img';
            spm_write_vol(VO,maskdat);
        end

        maskInfo = iimg_read_img('mask.img',2); % get info struct w/list of in-mask vox

        % define mat file names
        n = size(imagenames,1);
        if i_start==1
            for i = 1:n
                matnames{i} = fullfile(pwd,[subjnames{i} '_sn3d.mat']);
            end
        else
            for i = 1:n
                matnames{i} = fullfile([pwd filesep 'iteration' num2str(i_start-1)],[subjnames{i} '_sn3d.mat']);
            end
        end

        % map image volumes into memory
        disp('Mapping images to warp');
        VF = spm_vol(imagenames);

    end


    % ___________________________________________________________________
    %
    % Inline: Setup Normalization options
    % ___________________________________________________________________

    function setup_spm_norm_options(dographics)
        flags = struct('smosrc',0,'smoref',0,'regtype','mni',...
            'cutoff',30,'nits',16,'reg',0.1,'graphics',dographics);
        wflags=struct('wrap',[0 0 0],'vox',Inf,'bb',Inf,'preserve',0);
        % VG        - template handle(s)
        % VF        - handle of image to estimate params from
        % matname   - name of file to store deformation definitions
        % VWG       - template weighting image
        % VWF       - source weighting image
        % flags     - flags.  If any field is not passed, then defaults are assumed.
        %             smosrc - smoothing of source image (FWHM of Gaussian in mm).
        %                      Defaults to 8.
        %             smoref - smoothing of template image (defaults to 0).
        %             regtype - regularisation type for affine registration
        %                       See spm_affreg.m (default = 'mni').
        %             cutoff  - Cutoff of the DCT bases.  Lower values mean more
        %                       basis functions are used (default = 30mm).
        %             nits    - number of nonlinear iterations (default=16).
        %             reg     - amount of regularisation (default=0.1)

        % weighting for source smoothing
        % now set to: 16*i^-1.5
        wfun = inline('i .^ -a','i','a');
        source_w =max([16.*wfun(1:iter,0.7); ones(1,iter)*2],[],1);

        cutoff_w = 30 .* wfun(1:iter,.3);

        fprintf(1,'Initial source smoothing: %3.2f mm.  Source smoothing after %3.0f iterations: %3.2f mm\n',source_w(1),iter,source_w(end));
        fprintf(1,'Initial basis function smallest period: %3.2f mm.  After %3.0f iterations: %3.2f mm\n',cutoff_w(1),iter,cutoff_w(end));

    end




    % ___________________________________________________________________
    %
    % Inline: Get name of output param file and save eval string for 1
    % subject
    % ___________________________________________________________________
    function [VG, save_subj_str] = get_savename;
        VG = spm_vol('mean.img');
        [dd,ff] = fileparts(prmnames(i,:));
        save_subj_str = ['save ' ff ' -append VG Tr'];
        disp('Will save params for this subject to:');
        disp(save_subj_str)
    end



end %% end main function


