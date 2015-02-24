function EXPT = batch_run_rfx(EXPT)
% function EXPT = batch_run_rfx(EXPT)
% runs correlation analyses with behavior over whole brain
% using SPM
%
% by Tor Wager

BAT.estnow = 1;     % estimate the model
BAT.estcons = 1;    % estimate contrasts
BAT.dores = 1;      % get results

BAT.M_I = 0;

if ~isfield(EXPT,'behavior')
    EXPT.behavior = [];
end

if ~isfield(EXPT.RFX,'mask')
    EXPT.RFX.mask = [];
end

BAT.covariate = EXPT.behavior;



if ~isfield(EXPT.RFX,'connums')
    disp(['Making EXPT.RFX.connums from EXPT.SNPM'])
    EXPT.RFX.connums = EXPT.SNPM.connums;
end

if ~isfield(EXPT.RFX,'connames')
    disp(['Making EXPT.RFX.connames from EXPT.SNPM'])
    EXPT.RFX.connames = EXPT.SNPM.connames;
end
    
if ~isfield(EXPT.RFX,'mask') & isfield(EXPT.SNPM,'mask')
    disp(['Making EXPT.RFX.mask from EXPT.SNPM'])
    EXPT.RFX.mask = EXPT.SNPM.mask;
end

if ~isfield(EXPT.RFX,'P') & isfield(EXPT.SNPM,'P')
    disp(['Making EXPT.RFX.P from EXPT.SNPM'])
    EXPT.RFX.P = EXPT.SNPM.P;
end

BAT.multCompCorrect = EXPT.RFX.multCompCorrect; 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

BAT.u = EXPT.RFX.u;
	% height threshold - T or p value

BAT.k = EXPT.RFX.k;
	% extent threshold - number of voxels

interact = input('Interactive selection? (1/0) ');
wh = EXPT.RFX.connums;
if interact, BAT.pause = 1;, else, BAT.pause = 0;, wh=input('Enter which cons to run: ');, end
    
    
conindex = 1;

for mycon = wh

	conindex  = find(EXPT.RFX.connums == mycon);

    if mycon < 10, myz = '000';, else myz = '00';, end
    
    if ~(exist(['rfx' myz num2str(mycon)],'dir'))
            eval(['mkdir rfx' myz num2str(mycon)])
    end
    
    eval(['cd rfx' myz num2str(mycon)])
    
    %EXPT.RFX.mask = EXPT.SNPM.mask; % {['../rfx' myz num2str(mycon) '/res_mask.img']};
    %explicit mask to use--goes to M_P
    %type in path
    % e.g., TOR.MaskImg = {'mask.img'};
    % use {} if no mask to be used
    
    if size(EXPT.RFX.mask,1) > 1
        BAT.MaskImg = {[deblank(EXPT.RFX.mask(conindex,:))]};
        BAT.expMask = 1;
    elseif isfield(EXPT.RFX,'mask')
        if ~isempty(EXPT.RFX.mask)
            BAT.MaskImg = {EXPT.RFX.mask};
            BAT.expMask = 1;
        else
            BAT.MaskImg = {''};
            BAT.expMask = 0;
        end
    else
        BAT.MaskImg = {''};
        BAT.expMask = 0;
    end
    
    BAT.P = EXPT.RFX.P{conindex};
    BAT.npairs = size(EXPT.RFX.P{conindex},1);
    
    a= deblank(EXPT.RFX.connames(conindex,:)); a(a == ' ') = '_';
    BAT.covname = ['corr_' a];  % myz num2str(mycon)];
    BAT.resdir = pwd;  

    pwd
    BAT.P
    BAT

    disp(['Testing contrast: ' a])
    
    if interact
        OK = input('Go OK? Press 1 to go or 0 to skip ');
    else
        OK = 1;
    end
    
    if OK, 
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
    

    if interact, 
            disp('Type return when ready to go on!')
        keyboard
        input('Press enter to continue...'),
    end
    
end