function EXPT = mask_first_subject(EXPT,applynow)
% function EXPT = mask_first_subject(EXPT,applynow)
% Tor Wager

if ~(exist('analysis_masks') == 7), mkdir analysis_masks, end

for conindex = 1:length(EXPT.SNPM.connums)

        if size(EXPT.SNPM.mask,1) > 1
            P = str2mat(EXPT.SNPM.P{conindex}(1,:),EXPT.SNPM.mask(conindex,:)); % mask should be 2nd image
        else
            P = str2mat(EXPT.SNPM.P{conindex}(1,:),EXPT.SNPM.mask);
        end
        
        
        % -------------------------------------------------------------------------
        % reslice the contrast img for first subject and replace name in EXPT
        % -------------------------------------------------------------------------
            [d,f,e] = fileparts(EXPT.SNPM.P{conindex}(1,:));
            if ~isfield(EXPT,'studydir'), EXPT.studydir = spm_get(-1,'*','Choose main analysis dir');, end
            Q = fullfile(EXPT.studydir,'analysis_masks',[EXPT.subjects{1} '_masked_' f e]);

            warning off
            if applynow,
                disp('Calculating mask using these images:'), P
                disp(['Masked contrast image written as ' Q])
                Q = spm_imcalc_ui(P,Q,'i1 .* (i2 > 0) + 0 ./ (i2 > 0)',{0 1 16 0});
            end
            warning on
            
            EXPT.SNPM.P{conindex} = str2mat(Q,EXPT.SNPM.P{conindex}(2:end,:));
            
        end
            
end
    

%EXPT.SNPM.mask = newMask;


return
