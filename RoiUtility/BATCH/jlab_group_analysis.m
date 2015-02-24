% jlab_group_analysis
% example script


d = spm_get(-1,'*','Choose main analysis dir');
d = fileparts(d);
eval(['cd(''' d ''')'])

disp(['jlab_group_analysis  <...> batch stats using spm'])
disp('_________________________________________________')
disp('Choose analysis type:')
disp('  1   -   SnPM')
disp('  2   -   correlation with behavior')
disp('  3   -   random effects model')
disp('  4   -   nonlinear fit SNPM rfx')
disp('  5   -   nonlinear fit random effects')
yourchoice = input('Which would you like?');


if ~(exist('EXPT.mat'))
    disp(['No EXPT file found: Creating.'])
    domask = input('Explicit mask for analysis? 0 = no, 1 = single mask, 2 = indiv mask for each contrast ');
    EXPT = get_expt_info(10,domask);
else
    load EXPT
end




switch yourchoice
    
case 1
    % SnPM analysis

    EXPT.SNPM.u = input('Enter the primary corrected alpha value: (e.g.,.05) ');
    EXPT.SNPM.assess_extent = input('Assess spatial extent? (1/0) ');
    EXPT.SNPM.extent_u = input('Enter primary threshold for cluster-level analysis(e.g., .01): ');

    save EXPT EXPT
    batch_run_snpm(EXPT)
    batch_snpm2clusters


case 2
    % simple correlation/regression
    
    EXPT.CORR.multCompCorrect = input('Enter multiple comparisons procedure choices are FWE|FDR|uncorrected ','s'); 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

    EXPT.CORR.u = input('Enter height threshold T or p (e.g., .05) ');
	% height threshold - T or p value

    EXPT.CORR.k = input('Enter cluster extent threshold (e.g., 20) ');
	% extent threshold - number of voxels
    
    save EXPT EXPT
    EXPT = batch_run_correl(EXPT);
  
    
case 3
    % random effects - one sample T-test in SPM
    
    EXPT.RFX.multCompCorrect = input('Enter multiple comparisons procedure choices are FWE|FDR|uncorrected ','s'); 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

    EXPT.RFX.u = input('Enter height threshold T or p (e.g., .05) ');
	% height threshold - T or p value

    EXPT.RFX.k = input('Enter cluster extent threshold (e.g., 20) ');
	% extent threshold - number of voxels
    
    save EXPT EXPT
    EXPT = batch_run_rfx(EXPT);
    
case 4
    % random effects analysis over nonlinear height, delay, etc. parameters using SnPM
    EXPT.SNPM.u = input('Enter the primary corrected alpha value: (e.g.,.05) ');
    EXPT.SNPM.assess_extent = input('Assess spatial extent? (1/0) ');
    EXPT.SNPM.extent_u = input('Enter primary threshold for cluster-level analysis(e.g., .01): ');

    save EXPT EXPT
    batch_run_nl_snpm(EXPT)
    batch_snpm2clusters
    
    
case 5
    % random effects analysis over nonlinear h,d,intercept
    EXPT.RFX.multCompCorrect = input('Enter multiple comparisons procedure choices are FWE|FDR|uncorrected ','s'); 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

    EXPT.RFX.u = input('Enter height threshold T or p (e.g., .05) ');
	% height threshold - T or p value

    EXPT.RFX.k = input('Enter cluster extent threshold (e.g., 20) ');
	% extent threshold - number of voxels
    
    save EXPT EXPT
    EXPT = batch_run_nl_rfx(EXPT);
    
otherwise
    error('Unknown analysis choice.')
end