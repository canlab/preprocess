function IN = jlab_stats(IN)
%
% IN is a structure with a number of fields.


% get subject code if it doesn't exist
% ====================================================================== 
if ~(isfield(IN,'SubjCode') == 1)
	IN.SubjCode = input('Enter subject code ','s');
end

if ~(isfield(IN,'taskcode') == 1)
	IN.taskcode = input('Enter task dir name or return for empty ','s');
end


% make image directory list and results directory
% ====================================================================== 
for i = 1:IN.nruns
    IN.imgdir{i} = fullfile(IN.studydir,IN.SubjCode,IN.taskcode,['scan' num2str(i)]);
end

IN.resdir = fullfile(IN.studydir,'RESULTS',IN.taskcode,IN.SubjCode);
if ~(exist([IN.studydir filesep 'RESULTS']) == 7), eval(['!mkdir ' IN.studydir filesep 'RESULTS']), end 
if ~(exist([IN.studydir filesep 'RESULTS' filesep IN.taskcode]) == 7), eval(['!mkdir ' IN.studydir filesep 'RESULTS' filesep IN.taskcode]), end 
if ~(exist(IN.resdir) == 7), eval(['!mkdir ' IN.resdir]), end 


% setup:
% Kalina Christoff's bigmask stuff for modifying analysis mask
% -------------------------------------------------------------------
% bigmask is created in jlab_preproc.m
IN.in_create = 'n';		% create bigmask image file
IN.in_segment = 'n';		% use already created segmented images
IN.in_modify = 'y';		% modify spmCFG.mat before estimating
IN.in_estim = 'y';			% estimate after modifying spmCFG.mat
        %IN.PG = OPT.canonicalT1;
	% template image (T1.img or T2.img)
    %IN.PF = fullfile(subjDir,'anatomy','nht1spgr.img');				
	% image to segment
IN.cfg_file = [IN.resdir filesep 'SPMcfg.mat'];
		% cfg.mat file to modify
IN.bigmaskimg = fullfile(IN.studydir,IN.SubjCode,'anatomy','bigmask.img');


% setup:
% other variables we need
% -------------------------------------------------------------------
IN.v = length(IN.c) ./ IN.nruns;		
	% number of conditions/trial types in each scan
 
IN.HParam = IN.HParam .* ones(1,IN.nruns);
    % create vector of Hparam values for each session
    
IN.numscans = IN.nruns;
    % for compatibility with tor_spm_fmri_spm_ui
    
% Setup: For estimating only a single contrast
% ======================================================================
if IN.SelectedContrast
	myContrastName = IN.con_names{IN.Ic};
	%IN.Ic = whichContrast;
else
    myContrastName = [];
end

% make sure we're in the correct directory
% ====================================================================== 
IN

eval(['cd ' IN.resdir])
cont = [];
if ~(isfield(IN,'goOK')), IN.goOK = 0; end
if ~(IN.goOK) == 1
	while ~(strcmp(cont,'y') | strcmp('cont','n'))
		cont = input(['Data will be saved in ' pwd '. Save results here (y\n)?'],'s');
		if strcmp(cont,'n'),error('Please change to the correct results directory.'),end
	end
end


% get the list of actual file names
% ====================================================================== 

for i = 1:length(IN.imgdir)
	[fNames,dummy] = spm_list_files(IN.imgdir{i},IN.imgwildcard);
	% ...and add the directory name
 	a = repmat([IN.imgdir{i} filesep],length(fNames),1);
	fNames = [a fNames];
	IN.imgnames{i} = fNames;
	if isempty(fNames),error(['No files matching ' IN.imgwildcard ' in ' IN.imgdir{i}]),end 
end
 
% save parameter structure as .mat file
% ====================================================================== 
save Stats_Input_Params IN   

   

% set up bigmask insertion stuff
% ====================================================================== 
if IN.dobigmask
	IN.bEstNow = 0;
end



% build and estimate the model (or just set up)
% ====================================================================== 
if IN.buildmodel
	tor_spm_fmri_spm_ui(IN)
end  


% do CFG modification for bigmask, and estimate if specified
% ====================================================================== 
if IN.dobigmask
	tor_glm_specmask(IN)
end



% get the contrasts and save in xCon.mat  
% ==============================================================================
if IN.estimatecontrasts
	estimateContrasts(IN.contrasts,IN.con_names,IN.con_type,IN.SelectedContrast,myContrastName)
end


% threshold and display results
% ====================================================================== 

% top part is for getting results of ALL contrasts
% bottom part (last line) is for getting results and choosing a contrast w/ the GUI
%	or choosing the single contrast number specified in IN.Ic

if IN.getresults

if isfield(IN,'Ic')
	if strcmp(IN.Ic,'all')
		load xCon
		for i = 1:length(xCon)
            % results in SPM window
            % -------------------------------------------------------
			IN.Ic = i;
			[hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(IN);
			spm_list('List',SPM,VOL,[],[],'',hReg);				% list volume stats in window
			spm_print
			pause(5)
            
            % results clusters, containing betas
            % -------------------------------------------------------
            if IN.getclusters
                [P,D] = spm_list_files(pwd,'beta*img');
                if ~isempty(SPM.XYZ) & any(IN.cluster_con == i)
                    clusters = tor_extract_rois(P,SPM,VOL);
                    eval(['save clust_con' num2str(i) ' clusters'])
                    
                    if IN.paramcheck
                        param_stability(IN,clusters)
                    end
                end
            end
		end

        % only works on unix/linux - 2nd line should be used for windows.
		try
            eval(['!mv spm99.ps ' IN.SubjCode '_results.ps'])
        catch
            eval(['!move spm99.ps ' IN.SubjCode '_results.ps'])
        end
	end
else

	[hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(IN);
	spm_list('List',SPM,VOL,[],[],'',hReg);							% list volume stats in window
end

end



return