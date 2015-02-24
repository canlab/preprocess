function EXPT = reslice_masks_and_update_expt(EXPT,varargin)
% EXPT = reslice_masks_and_update_expt(EXPT,[opt] reset 1st image name)
%
% Tor Wager
% run from main analysis directory
%
% EXPT is the structure created with get_expt_info
% A 2nd argument will change the name of the 1st subject's con*img file 
% for each contrast back to one comparable to the other subjects.
% Use this if you want to re-slice contrast masks.
disp('Mask updating for nonlinear fits: cond*height, delay, intercept.')
domasksmooth = input('Smooth masks?  Enter 0 for none, or FWHM in mm: ');
if ~domasksmooth, reslmask = input('Reslice masks now (1) or already in space of functionals (0)? ');,end
reslice1sts = input('Mask 1st subject and/or replace mask names? (1/0): ');
if reslice1sts, applynow = input('Masks already calculated (0) or need to be calculated now (1)? ');, end

EXPT.SNPM.domasksmooth = domasksmooth;

% -------------------------------------------------------------------------
%re-set mask names if necessary - just copy and paste into workspace
% -------------------------------------------------------------------------
doreset = 0; if length(varargin) > 0, doreset = 1;,end
if doreset
    for i = 1:length(EXPT.DX.nlconP)    % loop through height, delay, intercept
    disp(['Changing 1st image name in EXPT.DX.nlconP{i}{:}!'])
    for conindex = 1:length(EXPT.SNPM.connums)
        myfile = deblank(EXPT.DX.nlconP{i}{conindex}(2,:)); 
        [myd] = fileparts(deblank(EXPT.DX.nlconP{i}{conindex}(1,:))); 
        %[myd,f,e] = fileparts(fileparts(deblank(EXPT.SNPM.P{conindex}(1,:))));
        EXPT.DX.nlconP{i}{conindex} = str2mat([myd filesep myfile(end-11:end)],EXPT.SNPM.P{conindex}(2:end,:));
    end
    
    if ~(exist(EXPT.DX.nlconP{i}{1}(1,:)))
        myd = spm_get(-1,'con*img','Select dir containing 1st subject cond*height.img');
        for conindex = 1:length(EXPT.SNPM.connums)
            myfile = deblank(EXPT.DX.nlconP{i}{conindex}(2,:)); 
            EXPT.SNPM.P{conindex} = str2mat([myd filesep myfile(end-11:end)],EXPT.DX.nlconP{i}{conindex}(2:end,:));
        end
    end
    end %for i = 1:length(EXPT.DX.nlconP)    % loop through height, delay, intercept
end

if ~(exist('analysis_masks') == 7), mkdir analysis_masks, end

for i = 1:length(EXPT.DX.nlconP)    % loop through height, delay, intercept
    
    for conindex = 1:length(EXPT.SNPM.connums)

        if size(EXPT.SNPM.mask,1) > 1
            P = str2mat(EXPT.DX.nlconP{i}{conindex}(1,:),EXPT.SNPM.mask(conindex,:)); % mask should be 2nd image
        else
            P = str2mat(EXPT.DX.nlconP{i}{conindex}(1,:),EXPT.SNPM.mask);
        end
        
        
        % -------------------------------------------------------------------------
        % smooth the mask, if specified, and replace the mask name in EXPT.SNPM.mask
        % -------------------------------------------------------------------------

        if domasksmooth
           
            [d f e] = fileparts(P(2,:));
            
            Q = ['analysis_masks' filesep f '_smooth' num2str(domasksmooth) e];
            spm_smooth(P(2,:),Q,domasksmooth);
            P = str2mat(P(1,:),Q);
            
            if conindex == 1, newMask(1,:) = Q;
            else newMask = str2mat(newMask,Q);
            end
        else
            % -------------------------------------------------------------------------
            % reslice the mask and store in analysis_masks
            % -------------------------------------------------------------------------
            if reslmask
                nms = reslice_imgs(P(1,:),P(2,:),1);
                [d f e] = fileparts(P(2,:));
                outnm = fullfile(d,['r' f e]); 
                if conindex == 1, newMask(1,:) = outnm;
                else newMask = str2mat(newMask,outnm);
                end
            else
                outnm = P(2,:);
                if conindex == 1, newMask(1,:) = outnm;
                else newMask = str2mat(newMask,outnm);
                end
            end   
        end
        
        % -------------------------------------------------------------------------
        % reslice the contrast img for first subject and replace name in EXPT
        % -------------------------------------------------------------------------
        if reslice1sts
            [d,f,e] = fileparts(EXPT.DX.nlconP{i}{conindex}(1,:));
            Q = fullfile(EXPT.studydir,'analysis_masks',[EXPT.subjects{1} '_masked_' f e]);
            disp('Calculating mask using these images:'), P
            disp(['Masked contrast image written as ' Q])
            warning off
            if applynow,
                Q = spm_imcalc_ui(P,Q,'i1 .* (i2 > 0) + 0 ./ (i2 > 0)',{0 1 16 0});
            end
            EXPT.DX.nlconP{i}{conindex} = str2mat(Q,EXPT.DX.nlconP{i}{conindex}(2:end,:));
            warning on
        end
        
    end
end
    

EXPT.SNPM.mask = newMask;