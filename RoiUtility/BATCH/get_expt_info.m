function EXPT = get_expt_info(varargin)
% EXPT = get_expt_info([opt] EXPT)
%
% Start running in the directory above individual results directories!!
%
% compiles info about your experiment for use with other scripts
% e.g., batch_multisubject_cluster_ui
% batch_snpm2clusters
% excludes directories called rfx*, fix*, or snpm*, as these are likely to
% contain random effects results.
%
% also compiles info for random effects spm/snpm analysis, stored in
% EXPT.SNPM.  P{i} contains image file names for each contrast
% where i indexes contrasts.  Each contrast is an rfx analysis.
%
% dxTRs is the number of time points to estimate in deconvolution matrix
% snpm_explicit_mask: 2, 1 or 0, to select a mask image and apply it to the
% first subject, for restricted volume snpm analysis.  2 is individual masks
% for each contrast.  1 is same mask for all contrasts.
%
% Tor Wager, 04/04/02

clear EXPT

if length(varargin) > 0, 
    EXPT = varargin{1};
end

EXPT.studydir = spm_get(-1,'*','Select main analysis dir (contains ind sub analyses)');
eval(['cd ' EXPT.studydir])
disp(['Main study directorym (EXPT.studydir) is: ' EXPT.studydir]);

if ~isfield(EXPT,'subdir'), EXPT.subjects = [];,end
if isempty(EXPT.subjects)
    wc = input('No subject codes (EXPT.subjects) found.  Enter wildcard (e.g., 04*) :','s');
    d = dir(wc);
    for i = 1:length(d), EXPT.subjects{i} = d(i).name;,end
end
disp(['Subjects:'])
EXPT.subjects

% -------------------------------------------------------------------------------------------------
% * MENU
% -------------------------------------------------------------------------------------------------

fprintf(1,'Menu selections\n')
fprintf(1,'0\tGet names of anatomical and functional raw images for each subject (for ind. stats)\n')
fprintf(1,'1\tGet image and analysis files from individual sub results directories\n')
fprintf(1,'2\tSelect mask image(s) - either a single one or one for each contrast\n')
fprintf(1,'3\tSmooth mask image(s)\n')
fprintf(1,'4\tReslice mask image(s) to be in space of functionals\n')
fprintf(1,'5\tApply mask to first subject and create new masked con*img files AND replace names\n')
fprintf(1,'6\tReplace con*img names in experiment with masked con*img names (replace names ONLY)\n')
fprintf(1,'7\tBuild deconvolution design and add to EXPT\n')
fprintf(1,'8\tApply masks, etc. to deconv design\n')
fprintf(1,'9\tGet img and analysis files for existing EXPT.subjects in order\n')
mychoice = input('Enter vector of numbers for operations you want to perform, e.g., [1 3 4]: ');

if any(mychoice == 0)
    EXPT = getfunctnames(EXPT);
end

if any(mychoice == 1)
    EXPT = read_ind_results(EXPT);
end

 if any(mychoice == 9)
    EXPT = read_ind_results_existsubdir(EXPT);
 end

 

if any(mychoice == 2)
    snpm_explicit_mask = input('Press 1 for single mask, or 2 to choose masks for each contrast: ');
    EXPT = enter_EXPT_masks(EXPT,snpm_explicit_mask);
end

if any(mychoice == 3)
    EXPT = smooth_masks(EXPT);    
end
  
if any(mychoice == 4)
    EXPT = reslice_masks(EXPT);    
end

if any(mychoice == 5)
    EXPT = mask_first_subject(EXPT,1);    
end

if any(mychoice == 6)
    EXPT = mask_first_subject(EXPT,0);   
end

 if any(mychoice == 8)
    if doreslice & ~isempty(EXPT.DX.nlconP), EXPT = reslice_nl_masks_and_upd_expt(EXPT);,end
 end
 

 
if any(mychoice == 7)
    EXPT = tor_build_deconv_design_ui(EXPT);,
end

    % options for running single subject stats
if isfield(EXPT,'SUBJECT')
     
    if ~isfield(EXPT.SUBJECT,'studydir'),
        EXPT.SUBJECT.studydir = spm_get(-1,'*','Choose main expt dir');
    end
    if ~isfield(EXPT.SUBJECT,'object'),
        EXPT.SUBJECT.object = spm_get(1,'*img','Choose normalization object (T1spgr)');
    end
    if ~isfield(EXPT.SUBJECT,'canonicalT1'),
        EXPT.SUBJECT.canonicalT1 = spm_get(1,'*img','Choose normalization template');
    end
 
end
 %save EXPT EXPT
 
 
 
 return