function EXPT = smooth_masks(EXPT)

for conindex = 1:length(EXPT.SNPM.connums)

        if size(EXPT.SNPM.mask,1) > 1
            P = str2mat(EXPT.SNPM.P{conindex}(1,:),EXPT.SNPM.mask(conindex,:)); % mask should be 2nd image
        else
            P = str2mat(EXPT.SNPM.P{conindex}(1,:),EXPT.SNPM.mask);
        end
        
        
        % -------------------------------------------------------------------------
        % smooth the mask, if specified, and replace the mask name in EXPT.SNPM.mask
        % -------------------------------------------------------------------------
           
            [d f e] = fileparts(P(2,:));
            
            Q = ['analysis_masks' filesep f '_smooth' num2str(domasksmooth) e];
            spm_smooth(P(2,:),Q,domasksmooth);
            P = str2mat(P(1,:),Q);
            
            if conindex == 1, newMask(1,:) = Q;
            else newMask = str2mat(newMask,Q);
            end
            
end
    

EXPT.SNPM.mask = newMask;


return
