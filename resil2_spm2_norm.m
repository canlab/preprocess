function EXPT = resil2_spm2_norm(EXPT,docoreg,donorm,writeanat,writefunc,dosmooth)
% EXPT = resil2_spm2_norm(EXPT,docoreg,donorm,writeanat,writefunc,dosmooth)
% 
% donorm:   1 or 0, do normalization of anatomical and write
% writenorm: write normalized functional images
% dosmooth: smooth all functional images

template_file = which('scalped_avg152T1.img');
template_handle = spm_vol(template_file);

EXPT.FILES.spgr = [];

spm_defaults

anatvoxsize = [3.125 3.125 3]; %[.94 .94 1.5];
funcvoxsize = [3.125 3.125 3];
fwhm = [8 8 8];                 % smoothing kernel

for s = 1:length(EXPT.subjects)
    
    disp(['Subject ' num2str(s) ' : ' EXPT.subjects{s}])
    
    p = [EXPT.subjects{s} filesep 'anatomy/het1spgr.img'];
    
    anat_handle = spm_vol(p);
    matname = [spm_str_manip(p,'sd') '_sn.mat'];

    EXPT.FILES.spgr{s} = p;
    
    if docoreg
        
        disp('  Getting first functional image for coregistration')
        for i=1, 
            pp=spm_get('Files',fullfile(EXPT.subjects{s},'func','resilience',['run_' num2str(i)],'ra_img'),'ra*0001.img');,
        end
    
        EXPT.FILES.firstfunct{s} = pp;
        
        scnlab_coreg_anat2funct(EXPT.FILES.spgr{s},pp);
        
    end
        
        
        
    if donorm
       % anatomical normalization
        % ---------------------------------------------
    
        % Determine parameters
        spm_normalise(template_handle, anat_handle, matname, '', '', defaults.normalise.estimate);
    end
    
    if writeanat
        % write anatomical
        defaults.normalise.write.vox = anatvoxsize;
        spm_write_sn(anat_handle, matname, defaults.normalise.write);
    end
    

    if writefunc
        % get file names for functional normalization (creates w* in spm2)
        % ---------------------------------------------
        pp=[]; 
     
        disp('  Getting functional files for normalization')
        for i=1:6, 
            pp=str2mat(pp,spm_get('Files',fullfile(EXPT.subjects{s},'func','resilience',['run_' num2str(i)],'ra_img'),'ra*.img'));,
        end

        pp(1,:) = [];
    
        EXPT.FILES.ra_imgs{s} = pp;
    
        % write functionals
        % ---------------------------------------------
        defaults.normalise.write.vox = funcvoxsize;
        spm_write_sn(pp, matname, defaults.normalise.write);
    end
    
    
    if dosmooth
     
    % get handles for interactive window
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Smooth');

    % get file names for smoothing (w* in spm2, n* in spm 99)
    % ---------------------------------------------
    pp=[]; allpp = [];
     
    disp('  Getting functional files for smoothing')
    for i=1:6,  % loop through runs
        
        
        pp = spm_get('Files',fullfile(EXPT.subjects{s},'func','resilience',['run_' num2str(i)],'ra_img'),'wra*.img');
        allpp=str2mat(allpp,pp);,
    

        % do the smoothing
        % ---------------------------------------------
        spm('Pointer','Watch');
        spm('FigName','Smooth: working',Finter,CmdLine);
        spm_progress_bar('Init',size(pp,1),'Smoothing','Volumes Complete');
        for j = 1:size(pp, 1)
            curr = deblank(pp(j,:));                                % input image
            [path,name,ext,ver] = fileparts(curr);                  % old path - same as input
            path = fullfile(EXPT.subjects{s},['scan' num2str(i)]);    % new path - scan directory
            smoothed_name = fullfile(path,['s' name ext ver]);      % output image
            spm_smooth(curr,smoothed_name,fwhm);
            spm_progress_bar('Set',j);
        end
        spm_progress_bar('Clear',j);
        spm('FigName','Smooth: done',Finter,CmdLine);
        spm('Pointer');

    end % loop through runs
    
    allpp(1,:) = [];
    EXPT.FILES.wra_imgs{s} = allpp;
          
    end % if dosmooth
        

end % subject loop

% FORMAT params = spm_normalise(VG,VF,matname,VWG,VWF,flags)
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
% ________________________________


% FORMAT VO = spm_write_sn(V,prm,flags,msk)
% V         - Images to transform (filenames or volume structure).
% matname   - Transformation information (filename or structure).
% flags     - flags structure, with fields...
%           interp   - interpolation method (0-7)
%           wrap     - wrap edges (e.g., [1 1 0] for 2D MRI sequences)
%           vox      - voxel sizes (3 element vector - in mm)
%                      Non-finite values mean use template vox.
%           bb       - bounding box (2x3 matrix - in mm)
%                      Non-finite values mean use template bb.
%           preserve - either 0 or 1.  A value of 1 will "modulate"
%                      the spatially normalised images so that total
%                      units are preserved, rather than just
%                      concentrations.
% msk       - An optional cell array for masking the spatially
%             normalised images (see below).
%
% Warped images are written prefixed by "w".
%
% Non-finite vox or bounding box suggests that values should be derived
% from the template image.
%
% Don't use interpolation methods greater than one for data containing
% NaNs.
