function EXPT = batch_run_correl(EXPT)
% function EXPT = batch_run_correl(EXPT)
% runs correlation analyses with behavior over whole brain
% using SPM
%
% by Tor Wager

BAT.estnow = 1;     % estimate the model
BAT.estcons = 1;    % estimate contrasts
BAT.dores = 1;      % get results

if ~isfield(EXPT,'behavior')
    EXPT.behavior = [];
end

if ~isfield(EXPT.RFX,'mask')
    EXPT.RFX.mask = [];
end

BAT.covariate = EXPT.behavior;

interact = input('Interactive selection? (1/0) ');
if interact, BAT.pause = 1;, else, BAT.pause = 0;, end

if ~isfield(EXPT.RFX,'connums')
    disp(['Making EXPT.RFX.connums from EXPT.SNPM'])
    EXPT.RFX.connums = EXPT.SNPM.connums;
end

if ~isfield(EXPT.RFX,'connames')
    disp(['Making EXPT.RFX.connames from EXPT.SNPM'])
    EXPT.RFX.connames = EXPT.SNPM.connames;
end
    
if ~isfield(EXPT.RFX,'mask') & isfield(EXPT.SNPM,'mask')
    disp(['Making EXPT.RFX.mask from EXPT.DX.nlconP'])
    EXPT.RFX.mask = EXPT.SNPM.mask;
end

%if ~isfield(EXPT.RFX,'P') & isfield(EXPT.DX,'nlconP')
%    disp(['Making EXPT.RFX.P from EXPT.DX.nlconP'])
%    EXPT.RFX.P = EXPT.DX.nlconP;
%end

BAT.multCompCorrect = EXPT.RFX.multCompCorrect; 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

BAT.u = EXPT.RFX.u;
	% height threshold - T or p value

BAT.k = EXPT.RFX.k;
	% extent threshold - number of voxels


    
    
conindex = 1;
overallconindex = 1;

for mycon = 1:size(EXPT.DX.contrasts,1)

    if mycon < 10, myz = '000';, else myz = '00';, end
    
    for myparam = 1:length(EXPT.DX.params)
        
            if ~(exist(['rfx' myz num2str(mycon) '_' EXPT.DX.params{myparam}],'dir'))
                eval(['mkdir rfx' myz num2str(mycon) '_' EXPT.DX.params{myparam}])
            end
    
            eval(['cd rfx' myz num2str(mycon) '_' EXPT.DX.params{myparam}])
         
    
    %EXPT.RFX.mask = EXPT.SNPM.mask; % {['../rfx' myz num2str(mycon) '/res_mask.img']};
    %explicit mask to use--goes to M_P
    %type in path
    % e.g., TOR.MaskImg = {'mask.img'};
    % use {} if no mask to be used
    
    if size(EXPT.RFX.mask,1) > 1
        BAT.MaskImg = {deblank(EXPT.RFX.mask(overallconindex,:))}; % {['../' deblank(EXPT.RFX.mask(conindex,:))]};
        BAT.expMask = 1;
    elseif isfield(EXPT.RFX,'mask')
        if ~isempty(EXPT.RFX.mask)
            BAT.MaskImg = {[EXPT.RFX.mask]}; % {['../' EXPT.RFX.mask]};
            BAT.expMask = 1;
        else
            BAT.MaskImg = {''};
            BAT.expMask = 0;
        end
    else
        BAT.MaskImg = {''};
        BAT.expMask = 0;
    end
    
    BAT.P = EXPT.DX.nlconP{myparam}{mycon};
    BAT.npairs = size(BAT.P,1);
    
    delayOK = 1;
    
    % for testing delay only in regions with significant height
    % --------------------------------------------------------------
    %if myparam == 2     % delay parameter
        % special: if delay (param 2), only run if height was significant!  
        % and mask with height significant mask
    %    heightmask = [BAT.resdir filesep 'res_mask.img'];
    %    BAT.MaskImg = heightmask;
    %    BAT.expMask = 1;
    %    V = spm_vol(heightmask);
    %    vol = spm_read_vols(V);
    %    if ~any(any(any(vol))), delayOK = 0;, end
    %end
    
    
    BAT.resdir = pwd;  

    pwd
    BAT.P
    BAT

    disp(['Testing contrast: ' num2str(mycon)])
    
    if interact
        OK = input('Go OK? Press 1 to go or 0 to skip ');
    else
        OK = 1;
    end
    
    if OK & delayOK, 
        tor_spm_random_effects(BAT), 
        % ----------------------------------------------------------------------
        % * EXTRACT CLUSTERS
        % ----------------------------------------------------------------------

        a = pwd; a = str2num(a(end-3:end));
        P = EXPT.RFX.P{EXPT.RFX.connums == a};
        [CLU,clusters] = spm2struct(pwd, EXPT.RFX, P); 

        if ~isempty(clusters), 
            save CLU clusters
            disp(['directory is ' pwd])
            conname = deblank(EXPT.RFX.connames(EXPT.RFX.connums == a,:));
            fprintf(1,'Contrast: %s',conname)
            cluster_table(clusters)
    
            for j = 1:length(clusters), clusters(j).title = conname;,end
            a = pwd; a = a(end-6:end);
            eval(['save ' a '_clusters clusters CLU']) 
        end

        % ----------------------------------------------------------------------    
    end
    
    cd ..
    conindex = conindex + 1;
    if interact, input('Press enter to continue...'),end
    
end % loop through params
overallconindex = overallconindex + 1;

end % loop through contrasts

return
