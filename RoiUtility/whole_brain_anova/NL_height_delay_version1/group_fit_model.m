function EXPT = group_fit_model(EXPT)
% function EXPT = group_fit_model(EXPT)
% Tor Wager 6/18/02

S = EXPT.FILT.S;

EXPT.DX.doridge = input('Use ridge regression (k=.0001)? enter 1/0 ');

% ------------------------------------------------------------------
% * set up mask, if not already done
% ------------------------------------------------------------------

Pm = input('Mask images and fit to only in-mask voxels? (1/0) ');

calcmask = 1;
if isfield(EXPT.DX,'mask')
    if ~isempty(EXPT.DX.mask)
        calcmask = 0;
    end
end

Pf = EXPT.im_files{1}(1,:);     % first functional image
if Pm & calcmask, 
    if ~isempty(EXPT.SNPM.mask), Pm = EXPT.SNPM.mask(1,:);, 
    else, Pm = spm_get(1,'*img','Select mask image');, 
    end
    fprintf(1,'\nMask image is %s\n',Pm)
    mask = get_mask_vol(Pf,Pm);
    EXPT.DX.mask = mask;
elseif ~(Pm), 
    EXPT.DX.mask = [];
end

% ------------------------------------------------------------------
% * extract betas and save image files
% ------------------------------------------------------------------  

for i = 1:length(EXPT.snums)
    
    if length(EXPT.DX.DX) > 1, 
        DX = EXPT.DX.DX{i};
    else
        DX = EXPT.DX.DX{1};
    end
    
    P = EXPT.im_files{i};
    
    if EXPT.runsubs(i)
        eval(['cd ' EXPT.subjects{i}])
    
        EXPT.dx_beta_imgs{i} = ind_model_fit(P,S,DX,EXPT.nsess,2,EXPT.DX.mask,EXPT.DX.doridge);
    
        cd ..
    end
     
 end
 