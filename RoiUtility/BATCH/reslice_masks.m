function EXPT = reslice_masks(EXPT)

newMask = [];

for conindex = 1:length(EXPT.SNPM.connums)

%for mymask =unique(EXPT.SNPM.mask,'rows')

        if size(EXPT.SNPM.mask,1) > 1
            P = str2mat(EXPT.SNPM.P{conindex}(2,:),EXPT.SNPM.mask(conindex,:)); % mask should be 2nd image
        else
            P = str2mat(EXPT.SNPM.P{conindex}(2,:),EXPT.SNPM.mask);
        end
        
	%P = str2mat(EXPT.SNPM.P{1}(1,:),mymask);
        
        % -------------------------------------------------------------------------
        % reslice the mask and store in analysis_masks
        % -------------------------------------------------------------------------
  		
                [d f e] = fileparts(P(2,:));
                outnm = fullfile(d,['r' f e]); 
		
		if isempty(strmatch(outnm,newMask))  
			disp(['RESLICING mask: ' P(2,:)])         
			nms = reslice_imgs(P(1,:),P(2,:),1);
		end


                if conindex == 1, newMask(1,:) = outnm;
                else newMask = str2mat(newMask,outnm);
                end
            
end
    

EXPT.SNPM.mask = newMask;


return
