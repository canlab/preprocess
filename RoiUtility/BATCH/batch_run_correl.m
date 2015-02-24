function EXPT = batch_run_correl(EXPT)
% function EXPT = batch_run_correl(EXPT)
% runs correlation analyses with behavior over whole brain
% using SPM
%
% by Tor Wager

BAT.estnow = 1;     % estimate the model
BAT.estcons = 1;    % estimate contrasts
BAT.dores = 1;      % get results

BAT.covariate = EXPT.behavior';
[mdir] = pwd;  
      
    
    % this option is for automatically taking significant activations for the contrast
    % as a mask for correlation analyses (reduces search space!)
    % 2 == use masks from nonlin fitting dirs, 1 == use masks from SNPM directories
use_snpm_res_masks = input('Use DX/SNPM/default analysis results masks for correlations? (2/1/0) ');


%if ~isfield(EXPT.CORR,'connums')
if use_snpm_res_masks == 1
    disp(['Making EXPT.CORR contrast info (connums, mask, connames) from EXPT.SNPM'])
    EXPT.CORR.connums = EXPT.SNPM.connums;
    EXPT.CORR.connames = EXPT.SNPM.connames;
    EXPT.CORR.P = EXPT.SNPM.P;
    %EXPT.CORR.mask = EXPT.SNPM.mask;
elseif use_snpm_res_masks == 2
    disp(['Making EXPT.CORR contrast info (connums, mask, connames) from EXPT.DX'])
    EXPT.CORR.connums = 1:size(EXPT.DX.contrasts,1);
    EXPT.CORR.connames = str2mat(EXPT.DX.connames);
    EXPT.CORR.P = EXPT.DX.nlconP{1};
else
    disp('Replacing EXPT.CORR.mask with EXPT.SNPM.mask contents.')
    EXPT.CORR.mask = EXPT.SNPM.mask;
end

%if ~isfield(EXPT.CORR,'connames')
%    disp(['Making EXPT.CORR.connames from EXPT.SNPM'])
%    EXPT.CORR.connames = EXPT.SNPM.connames;
%end
    
%if ~isfield(EXPT.CORR,'mask') & isfield(EXPT.SNPM,'mask')
%    disp(['Making EXPT.CORR.mask from EXPT.SNPM'])
%    EXPT.CORR.mask = EXPT.SNPM.mask;
%end

%if ~isfield(EXPT.CORR,'P') & isfield(EXPT.SNPM,'P')
%    disp(['Making EXPT.CORR.P from EXPT.SNPM'])
%    EXPT.CORR.P = EXPT.SNPM.P;
%end

BAT.multCompCorrect = EXPT.CORR.multCompCorrect; 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

BAT.u = EXPT.CORR.u;
	% height threshold - T or p value

BAT.k = EXPT.CORR.k;
	% extent threshold - number of voxels


BAT.M_I = 0;   % implicit masking - 1 or 0
                % I think implicit masking negates explicit masking
    

interact = input('Interactive selection? (1/0) ');
wh = EXPT.CORR.connums;
if interact, BAT.pause = 1;, else, BAT.pause = 0;, wh=input('Enter which cons to run: ');, end

conindex = 1;

for mycon = wh

	conindex  = find(EXPT.CORR.connums == mycon);

    if mycon < 10, myz = '000';, else myz = '00';, end
    
    if ~(exist(['corr' myz num2str(mycon)],'dir'))
            eval(['mkdir corr' myz num2str(mycon)])
    end
    
    eval(['cd corr' myz num2str(mycon)])
    
    %EXPT.CORR.mask = EXPT.SNPM.mask; % {['../rfx' myz num2str(mycon) '/res_mask.img']};
    %explicit mask to use--goes to M_P
    %type in path
    % e.g., TOR.MaskImg = {'mask.img'};
    % use {} if no mask to be used
    
    % this option is for automatically taking significant activations for the contrast
    % as a mask for correlation analyses (reduces search space!)
    % 2 == use masks from nonlin fitting dirs, 1 == use masks from SNPM directories
    if use_snpm_res_masks
        mymask = fullfile(mdir,['snpm' myz num2str(mycon)],'SnPMt_filtered.img');
        if use_snpm_res_masks == 2
            mymask = fullfile(mdir,['snpm' myz num2str(mycon) '_height'],'SnPMt_filtered.img');
        end
        
        if exist(mymask) == 2
            if conindex == 1
                EXPT.CORR.mask = mymask;
            else
                EXPT.CORR.mask = str2mat(EXPT.CORR.mask,mymask);
            end
        else
            mymask = EXPT.SNPM.mask(1,:);
            warning(['Mask not found: Using ' mymask])
          
            if conindex == 1
                EXPT.CORR.mask = mymask;
            else
                EXPT.CORR.mask = str2mat(EXPT.CORR.mask,mymask);
            end
        end
    end  
            
    
    if size(EXPT.CORR.mask,1) > 1
        BAT.MaskImg = {deblank(EXPT.CORR.mask(conindex,:))};
        BAT.expMask = 1;
    elseif isfield(EXPT.CORR,'mask')
        if ~isempty(EXPT.CORR.mask)
            BAT.MaskImg = {EXPT.CORR.mask};
            BAT.expMask = 1;
        else
            BAT.MaskImg = {''};
            BAT.expMask = 0;
        end
    else
        BAT.MaskImg = {''};
        BAT.expMask = 0;
    end
    
    BAT.P = EXPT.CORR.P{conindex};
    
    a= deblank(EXPT.CORR.connames(conindex,:)); a(a == ' ') = '_';
    BAT.covname = ['corr_' a];  % myz num2str(mycon)];
    BAT.resdir = pwd;  

    pwd
    BAT.P
    BAT

    
    if interact
        OK = input('Go OK? Press 1 to go or 0 to skip, or 2 to run all remaining');
        if OK == 2, interact = 0;, end
    else
        OK = 1;
    end
    
    if OK, tor_spm_simple_correl(BAT), end
    
    cd ..
    conindex = conindex + 1;
    if interact, input('Press enter to continue...'),end
    
end