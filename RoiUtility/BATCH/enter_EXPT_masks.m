function EXPT = enter_EXPT_masks(EXPT,snpm_explicit_mask)

EXPT.SNPM.mask = [];
if snpm_explicit_mask == 2, indivmasks = 1;, else, indivmasks = 0;,end

if ~indivmasks, 
	EXPT.SNPM.mask = spm_get(1,'*img','Select mask image.');,
	EXPT.SNPM.mask = repmat(EXPT.SNPM.mask,length(EXPT.SNPM.connums),1);

elseif indivmasks,

    for conindex = 1:length(EXPT.SNPM.connums)
       
            emask = spm_get(1,'*img',['Sel. mask for ' EXPT.SNPM.connames(conindex,:)]);,  
            if conindex == 1,
                EXPT.SNPM.mask(1,:) = emask;
            else
                EXPT.SNPM.mask = str2mat(EXPT.SNPM.mask,emask);
            end
    end

end



return
