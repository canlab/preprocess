function EXPT = group_brain_nlfit(EXPT)
% EXPT = group_brain_nlfit(EXPT)
%
% Start in main analysis directory, above individual subjects
%
% Tor Wager

S = EXPT.FILT.S;

[y,m,d,h,mi,s] = datevec(now);
str = (['group_brain_nlfit_' num2str(y) '_' num2str(m) '_' num2str(d) '.log'])
[FID, MESSAGE] = fopen(str,'a+');
disp(['Logged in ' str])
fprintf(FID,'Group_brain_nlfit > started %3.0f:%3.0f:%3.0f\n',h,mi,s)



% ------------------------------------------------------------------
% * set up input image list, if not already done
% ------------------------------------------------------------------
if ~isfield(EXPT,'dx_beta_imgs')
    index = 1;
    for i = EXPT.subjects
        eval(['cd ' i{1}])
        [EXPT.dx_beta_imgs{index},Dirs] = spm_list_files('.','dx_beta*img');
        if isempty(EXPT.dx_beta_imgs{index}), error(['No dx_beta imgs found in ' i{1}]), end
        index = index + 1;
        cd ..
    end
end



% ------------------------------------------------------------------
% * set up mask, if not already done
% ------------------------------------------------------------------

Pm = input('Mask images and fit to only in-mask voxels? (1/0) ');

calcmask = 1;
if isfield(EXPT.DX.mask)
    if ~isempty(EXPT.DX.mask)
        calcmask = 0;
    end
end

if exist(EXPT.dx_beta_imgs{1}(1,:))
    Pf = EXPT.dx_beta_imgs{1}(1,:); 
else    % try this
    Pf = [EXPT.subjects{1} filesep EXPT.dx_beta_imgs{1}(1,:)];     % first dx_beta image
    if ~exist(Pf)
	disp('dx_beta_imgs not found.  Please select any dx_beta_img.')
	Pf = spm_get(1,'dx_beta*img');
    end
end

if Pm & calcmask, 
    if ~isempty(EXPT.SNPM.mask), Pm = EXPT.SNPM.mask(1,:);, 
    else, Pm = spm_get(1,'*img','Select mask image');, 
    end
    fprintf(1,'\nMask image is %s\n',Pm)
    fprintf(FID,'\nMask image is %s\n',Pm)
    mask = get_mask_vol(Pf,Pm);
    EXPT.DX.mask = mask;
elseif ~(Pm), 
    EXPT.DX.mask = [];
end     

% ------------------------------------------------------------------
% * load betas and save height, delay, and intercept img files
% ------------------------------------------------------------------  

for i = 1:length(EXPT.snums)
    
    %if length(EXPT.DX.DX) > 1, 
    %    DX = EXPT.DX.DX{i};
    %else
    %    DX = EXPT.DX.DX{1};
    %end
    
    if EXPT.runsubs(i)
    	P = EXPT.dx_beta_imgs{i};
    
    	fprintf(1,'\nStarting subject %s',EXPT.subjects{i})
	fprintf(3,'\nStarting subject %s',EXPT.subjects{i})
    	eval(['cd ' EXPT.subjects{i}])
	disp(['First image is ' P(1,:)])
       fprintf(3,'First image is %s', P(1,:))
    	EXPT.dx_nlfit_imgs{i} = ind_model_nlfit(P,[],EXPT.DX,2,EXPT.DX.mask,FID);
    
    	cd ..
    end
     
 end
 
fclose(FID) 

return

