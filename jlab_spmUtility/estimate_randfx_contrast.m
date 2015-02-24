clear c

cnum = 1;
c{cnum} = [1];
names{cnum} = pwd;
type{cnum} = 'T';

% define results thresholding choices
% [run tor_spm_results_ui(TOR)]
% ====================================================================== 

TOR.maskWOtherContrasts = 0;
	% mask with other contrasts; 1 or 0

TOR.multCompCorrect = 'uncorrected'; 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

TOR.u = .05;
	% height threshold - T or p value

TOR.k = 010;
	% extent threshold - number of voxels

TOR.Ic = 2;
	% optional: which contrast to show in results (number)
	% or 'all'.  If 'all', uses all contrasts, in order.

TOR.resdir = pwd;


estimateContrasts
[hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(TOR);
spm_list('List',SPM,VOL,[],[],'',hReg);	

% write results mask 
mask = voxel2mask(SPM.XYZ',VOL.DIM');
V = spm_vol('mask.img');   
V.fname = 'res_mask.img';
V = spm_write_vol(V,mask);  

