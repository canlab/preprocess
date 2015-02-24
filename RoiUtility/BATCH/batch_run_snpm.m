function batch_run_snpm(EXPT)
% function batch_run_snpm(EXPT)
% by Tor Wager

dosetup = EXPT.SNPM.dosetup;    % set up snpm cfg
docalc = EXPT.SNPM.dosetup;     % calculate perms
dores = 1;      % get results

interact = 0; %input('Interactive selection? (1/0) ');
if interact, 
	BAT.pause = 1;, 
else, 
	BAT.pause = 0;, 
	runcons = EXPT.SNPM.runcons; %input('Enter vector of contrasts to run or 0 for all:');
	if runcons == 0, runcons == EXPT.SNPM.connums;, end

end

if ~isfield(EXPT,'behavior')
    EXPT.behavior = [];
end

if ~isfield(EXPT.SNPM,'mask')
    domask = input('Explicit mask analysis? 1/0 ');
    if domask, EXPT.SNPM.mask = spm_get(1,'*.img','Choose analysis mask');,end
end

conindex = 1;

for mycon = EXPT.SNPM.connums

    if any(runcons == mycon)
    if mycon < 10, myz = '000';, else myz = '00';, end
    
    if ~(exist(['snpm' myz num2str(mycon)],'dir')) 
            eval(['mkdir snpm' myz num2str(mycon)])
    end
    
    eval(['cd snpm' myz num2str(mycon)])

    
    % for no masking: snpm will not do explicit masking.
    % to do it, you have to have images for at least 1 subject
    % masked with nan's.
    % if masking is specified,
    % take the first image, mask it, and replace the old name with this one in the analysis
    % for brain masking: images pre-masked with brain and surrounded with NaNs.
    % this is done in get_expt_info
    
    BAT.P = EXPT.SNPM.P{conindex};
    
    pwd
    BAT.P(1,:)
    EXPT.SNPM.connums
    EXPT.SNPM.mask(1,:)
    

            
    if interact
        OK = input('Go OK? Press 1 to go, 0 to skip, or 2 to run all remaining ');
        if OK == 2, interact = 0;, BAT.pause = 0;,end
    else
        if any(runcons == mycon), OK = 1;, else, OK = 0;, end
    end
    
    % snpm_batch_script

BAT.DesType =EXPT.SNPM.DesType;
    % Type of model
    % 'SingleSub: 2 conditions, replications',...
	% 'SingleSub: Single covariate of interest',...
	% 'MultiSub: 1 conditions, 1 scan per subject',...
	% 'MultiSub: 2 conditions, replications - randomization exp.',...
	% 'MultiSub: 2 conditions, replications - permutation test',...
	% 'MultiGroup: 2 groups, 2 conditions, 1 scan per condition',...
	% 'MultiGroup: 2 groups, 1 scan per subject',...
	% 'User Specified PlugIn',...

if BAT.DesType == 2, BAT.behavior = EXPT.cov;, EXPT.SNPM.numcovs = 0;, BAT.Xblk = EXPT.SNPM.Xblk; end

BAT.vFWHM = 10;
    % variance smoothing full-width half-max
    
BAT.bVolm = 1;
	% work volumetrically?

BAT.bST = 1;
    % collect suprathreshold stats?  1 = yes, 0 = no
    


BAT.iGloNorm = 1;
    % global normalization
	% '<no global normalisation>',...				%-1
	% 'proportional scaling',...				    %-2
	% 'AnCova',...						            %-3
	% 'AnCova {subject-specific}',...				%-4
	% 'AnCova {study-specific}');				    %-5
    
BAT.iGXcalc = 3;
    % global calculation
    %'omit',...							            %-1
    %'user specified',...					        %-2
    %'mean voxel value (within per image fullmean/8 mask)');	%-3

BAT.iTHRESH = 1;
    % threshold masking
    % 1 = none, 2 = proportional, 3 = absolute
    
BAT.iGMsca = 2; 
    % Grand mean scaling
    % 1 = scaling, 2 = none



% defined in plug-in
% ==============================================================

BAT.ConfCov = EXPT.SNPM.numcovs;
    % number of confounding covariates.  0-6.

   if ~isfield(EXPT.SNPM,'numcovs')
       % no covariates, do nothing
   elseif isempty(EXPT.SNPM.numcovs) | EXPT.SNPM.numcovs == 0
       %still do nothing
   else
       disp(['Scaling and using column vectors in EXPT.cov as covariates'])
       BAT.covs = scale(EXPT.cov);
       
       if BAT.DesType == 2, BAT.covs = scale(EXPT.SNPM.covnoint);,end
       
       disp(['Num covariates found in EXPT.cov: ' num2str(size(BAT.covs,2))])
       BAT.ConfCov = size(BAT.covs,2);
   end

BAT.bAproxTst = 0;
    % use approximate test? (1/0)

if 2^(size(BAT.P,1)) > 5000
    BAT.bAproxTst = 1;
    BAT.nperms = 5000;
    % nperms: How many permutations to do.  10000 is a good upper limit.
end

% defined in biman5_run_snpm
%BAT.P = getfiles2('../../imgs/*con_0002.img');
    % file names
    

% for results (batch_spm_snpm_pp)
% ==============================================================   
BAT.bNeg = EXPT.SNPM.negeffects == 1;
disp(['COMPUTING POS CONTRASTS (0) or NEG (1): ' num2str(BAT.bNeg)])

    % positive or negative effects
    % 0 = pos, 1 = neg

BAT.writeStat = 1;
    % write out statistic image? (1/0)
    
BAT.alpha = EXPT.SNPM.u;
    % t or p threshold (e.g., .05)
    
BAT.bSpatEx = EXPT.SNPM.assess_extent;
    % assess spatial extent? (1/0)
    
%BAT.clExtThresh = tinv(.995,size(BAT.P,1)-1);
BAT.clExtThresh = tinv(1-EXPT.SNPM.extent_u,size(BAT.P,1)-1);
    % primary t threshold for spatial extent
    
%BAT.pause = 1;
    % pause to view permutation distribution
    
if OK
    if dosetup, batch_spm_snpm_ui(BAT), end
    if docalc, spm_snpm(pwd), end
    if dores, batch_spm_snpm_pp(pwd,BAT), spm_print,end
end

!rm SnPM_ST.mat

    cd ..
    
    end % if runcons
        
    conindex = conindex + 1;
    
end % loop through contrasts
