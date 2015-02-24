function spm_snpm_pp(CWD,BAT)
% SnPM post processing and results display
% FORMAT spm_snpm_pp(CWD)
%
% CWD -	Directory containing SnPM results files
%
% If CWD is not specified then user is prompted to locate results file SnPM.mat
%_______________________________________________________________________
%
% spm_snpm_pp is the PostProcessing function for the SnPM nonParametric
% statistical analysis. SnPM statistical analyses are split into three
% stages; Setup, Compute & Assess. This is the third stage.
% Nonparametric randomisation distributions are read in from MatLab
% *.mat files, with which the observed statistic image is assessed
% according to user defined parameters. It is the SnPM equivalent of
% the "Results" section of SPM, albeit with reduced features.
%
% Voxel level corrected p-values are computed from the permutation
% distribution of the maximal statistic. If suprathreshold cluster
% statistics were collected in the computation stage (and the large
% SnPM_ST.mat file hasn't been deleted!), then assessment by
% suprathreshold cluster size is also available, using a user-specified
% primary threshold.
%
% Instructions:
%=======================================================================
%
% You are prompted for the following:
%
% (1) ResultsDir: If the results directory wasn't specified on the command
%     line, you are prompted to locate the SnPM results file SnPM.mat.
%     The directory in which this file resides is taken to be the
%     results directory, which must contain *all* the files listed
%     below ("SnPM files required").
%
%     Results (spm.ps & any requested image files) are written in the
%     present working directory, *not* the directory containing the
%     results of the SnPM computations.
%
%                           ----------------
%
% (2) +/-: Having located and loaded the results files, you are asked to
%     chose between "Positive or negative effects?". SnPM, like SPM,
%     only implements single tailed tests. Choose "+ve" if you wish to
%     assess the statistic image for large values, indicating evidence
%     against the null hypothesis in favour of a positive alternative
%     (activation, or positive slope in a covariate analysis).
%
%     Choose "-ve" to assess the negative contrast, i.e. to look for
%     evidence against the null hypothesis in favour of a negative
%     alternative (de-activation, or a negative slope in a covariate
%     analysis). The "-ve" option negates the statistic image and
%     contrast, acting as if the negative of the actual contrast was
%     entered.
%
%     A two-sided test may be constructed by doing two separate
%     analyses, one for each tail, at half the chosen significance
%     level, doubling the resulting p-values.
%     ( Strictly speaking, this is not equivalent to a rigorous two-sided )
%     ( non-parametric test using the permutation distribution of the     )
%     ( absolute maximum statistic, but it'll do!                         )
%
%                           ----------------
%
% (3) WriteFiles: After a short pause while the statistic image SnPMt.mat
%     is read in and processed, you have the option of writing out the
%     complete statistic image as an Analyze file.
%
%     The image is written to the present working directory, as
%     SPMt.{hdr,img} or SPMt_neg.{hdr,img}, "_neg" being appended when
%     assessing the -ve contrast. (So SPMt_neg.{hdr,img} is basically
%     the inverse of SPMt.{hdr,img}!).   This, and all othet images created
%     have "float" precision.  All voxels surviving the "Grey Matter
%     threshold" are written, remaining image pixels being given zero value.
%
%     Next you are asked if you want to write a filtered image.  This
%     images is not written immediately, but after the processing.  If there
%     are any significant voxels, an image will be created consisting soley
%     of significant regions; all other voxels will be zero.  You are asked
%     to name this file; it will be written to the present working directory.
%
%     Next you are given the option of writing the complete
%     single-step adjusted p-value image. This image has voxel values
%     that are the non-parametric corrected p-value for that voxel (the
%     proportion of the permutation distribution for the maximal
%     statistic which exceeds the statistic image at that voxel). Since
%     a small p indicates strong evidence, 1-p is written out. So,
%     large values in the 1-p image indicate strong evidence against
%     the null hypothesis.  The image is written to the
%     present working directory, as SnPMp_SSadj.{img,hdr} or
%     SnPMp_SSadj_neg.{img,hdr}. This image is computed on the fly, so
%     there may be a slight delay...
%
%     Note that voxel level permutation distributions are not
%     collected, so "uncorrected" p-values cannot be obtained.
%
%                           ----------------
%
% Next come parameters for the assessment of the statistic image...
%
% (4) alpha: (Corrected p-value for filtering)
%     Next you enter the \alpha level, the statistical significance
%     level at which you wish to assess the evidence against the null
%     hypothesis. In SPM this is called "filtering by corrected
%     p-value".  SnPM will only show you voxels (& suprathreshold
%     regions if you choose) that are significant (accounting for
%     multiple comparisons) at level \alpha.  I.e. only voxels (&
%     regions) with corrected p-value less than \alpha are shown to
%     you.
%
%     Setting \alpha to 1 will show you all voxels with a positive statistic.
%
%
% (5) SpatEx: If you collected supra-threshold cluster statistics during
%     the SnPM computation phase, you are offered the option to assess
%     the statistic image by supra-threshold cluster size (spatial
%     extent).
%
% 5a) ST_Ut: If you chose to asses spatial extent, you are now prompted
%     for the primary threshold. This is the threshold applied to the
%     statistic image for the identification of supra-threshold
%     clusters.
%
%     The acceptable range is limited.  SnPM has to collect
%     suprathreshold information for every relabelling. Rather that
%     pre-specify the primary threshold, information is recorded for
%     each voxel exceeding a low threshold (set in spm_snpm) for every
%     permutation. From this, suprathreshold cluster statistics can be
%     generated for any threshold higher than the low recording
%     threshold. This presents a lower limit on the possible primary
%     threshold.
%
%     The upper limit (if specified) corresponds to the statistic value
%     at which voxels become individually significant at the chosen
%     level (\alpha).  There is little point perusing a suprathreshold
%     cluster analysis at a threshold at which the voxels are
%     individually significant.
%
%     If the statistics are t-statistics, then you can also specify the
%     threshold via the upper tail probability of the t-distribution.
%
%     (NB: For the moment, \alpha=1 precludes suprathreshold analysis, )
%     (    since all voxels are significant at \alpha=1.               )
%
%
% That's it. SnPM will now compute the appropriate significances,
% reporting its progress in the MatLab command window. Note that
% computing suprathreshold cluster size probabilities can take a long
% time, particularly for low thresholds or large numbers of
% relabellings. 
%
% Eventually, the Graphics window will come up and the results displayed.
% First the histograms of the permutation distributions are displayed; these
% distributions determined the corrected p-values and are central to each
% data analysis.  They represent the typical (maximal) values under the null
% hypothesis.  The correctly labeled or "observed" value is noted, as well
% as the threshold; if a corrected p-value threshold of 0.05 is used then
% the threhold falls at the 95%-ile of the permutation distribution.
%
% Press RETURN in the command window to see the standard results page.
%     
% - Results
%=======================================================================
%
% The format of the results page is similar to that of SPM:
%
% A Maximum Intensity Projection (MIP) of the statistic image is shown
% top left: Only voxels significant (corrected) at the chosen level
% \alpha are shown. (If suprathreshold cluster size is being assessed,
% then clusters are shown if they have significant size *or* if they
% contain voxels themselves significant at the voxel level.) The MIP is
% labelled SnPM{t} or SnPM{Pseudo-t}, the latter indicating that
% variance smoothing was carried out.
%
% On the top right a graphical representation of the Design matrix is
% shown, with the contrast illustrated above.
%
% The lower half of the output contains the table of p-values and
% statistics, and the footnote of analysis parameters. As with SPM, the
% MIP is tabulated by clusters of voxels, showing the maximum voxel
% within each cluster, along with at most three other local maxima
% within the cluster. The table has the following columns:
% 
% * region: The id number of the suprathreshold. Number 1 is assigned to
%   the cluster with the largest maxima.
%
% * size{k}: The size (in voxels) of the cluster.
%
% * P(Kmax>=k): P-value (corrected) for the suprathreshold cluster size.
%   This is the probability (conditional on the data) of the experiment
%   giving a suprathreshold cluster of size as or more extreme anywhere
%   in the statistic image. This is the proportion of the permutation
%   distribution of the maximal suprathreshold cluster size exceeding
%   (or equalling) the observed size of the current cluster. This
%   field is only shown when assessing "spatial extent".
%
% * t / Pseudo-t: The statistic value.
%
% * P(Tmax>=t): P-value (corrected) for the voxel statistic.
%   This is the probability of the experiment giving a voxel statistic
%   this extreme anywhere in the statistic image. This is the
%   proportion of the permutation distribution of the maximal
%   suprathreshold cluster size exceeding (or equalling) the observed
%   size of the current cluster.
%
% * (uncorrected): If the statistic is a t-statistic (i.e. variance smoothing
%   was *not* carried out, then this field is computed by comparing the
%   t-statistic against Students t-distribution of appropriate degrees
%   of freedom. Thus, this is a parametric uncorrected p-value. If
%   using Pseudo-t statistics, then this field is not shown, since the
%   data necessary for computing non-parametric uncorrected p-values is
%   not computed by SnPM.
%
% * {x,y,z} mm: Locations of local maxima.
%
% The SnPM parameters footnote contains the following information:
%
% * Primary threshold: If assessing "spatial extent", the primary
%   threshold used for identification of suprathreshold clusters is
%   printed. If using t-statistics (as opposed to Pseudo-t's), the
%   corresponding upper tail probability is also given.
%
% * Critical STCS: The critical suprathreshold cluster size. This is
%   size above which suprathreshold clusters have significant size at
%   level \alpha It is computed as the 100(1-alpha)%-ile of the
%   permutation distribution of the maximal suprathreshold cluster
%   size. Only shown when assessing "spatial extent".
%
% * alpha: The test level specified.
%
% * Critical threshold: The critical statistic level. This is the value 
%   above which voxels are significant (corrected) at level \alpha.  It
%   is computed as the 100(1-alpha)%-ile of the permutation
%   distribution of the maximal statistic.
%
% * df: The degrees of freedom of the t-statistic. This is printed even
%   if
%   variance smoothing is used, as a guide.
%
% * Volume & voxel dimensions:
%
% * Design: Description of the design
%
% * Perms: Description of the exchangability and permutations used.
%
%
%
% - SnPM files required:
%=======================================================================
% spm_snpm_pp loads parameters and results from the following files, which
% must all be in the same directory:
%       SnPMcfg.mat    - SnPM design configuration
%       SnPM.mat        - SnPM analysis & permutation distribution
%       SnPMt.mat       - Pointlist of (Pseudo) t-statisic for actual labelling
%       XYZ.mat         - Co-ordinates of pointlist
%       SnPM_ST.mat (*) - Suprathreshold cluster statistics (if required)
%
% Further details of the actual variables required from these files are
% given in the main body of spn_snpm_pp
%
% (*) The SnPM_ST.mat file containing the suprathreshold cluster
% information for each of the relabellings can be very large, and is
% only needed if a suprathreshold cluster size test is required. If
% such an analysis is not required, but suprathreshold cluster stats
% were collected, then this file may be deleted, without compromising
% further voxel-level analyses.
%
%_______________________________________________________________________
% @(#)spm_snpm_pp.m	2.5 Andrew Holmes 02/01/17

%-----------------------------functions-called------------------------
% spm_DesMtx
% spm_Tcdf
% spm_clf
% spm_clusters
% spm_figure
% spm_get
% spm_hwrite
% spm_input
% spm_invTcdf
% spm_max
% spm_mip
% spm_snpm_pp
% spm_str_manip
% spm_t2z
% spm_type
% spm_xyz2e
%-----------------------------functions-called------------------------

%-Variable "decoder" - Following files/variables are required:
%=======================================================================
% NB: Mat files contain additional variables beyond those required here.
%     See function that wrote each file for full definitions.
% SnPM design configuration file: 			SnPMcfg.mat
%-----------------------------------------------------------------------
% - saved by spm_snpm_ui
% H             condition partition of DesMtx for correctly labeled data
% C             covariate partition of DesMtx for correctly labeled data
% B             block     partition of DesMtx for correctly labeled data
% G             confound  partition of DesMtx for correctly labeled data
% HCBGnames     string matrix of column names of [H C B G]
% PiCond        Permuted conditions matrix, one labelling per row, actual
%               labelling on first row
% sPiCond       String describing permutations in PiCond
% bhPerms       Flag indicating use of "half permutations" trick
% CONT          Contrast (only one)
% bVarSm        Flag for variance smoothing (Pseudo t-statistics)
% sVarSm        Sring describing variance Smoothing (empty if bVarSm=0)
% bST           Flag for collection of superthreshold info 
% sDesign       Description of PlugIn design
% 
% SnPM analysis & permutation distribution file:	SnPM.mat
%-----------------------------------------------------------------------
% - saved by spm_snpm
% S		- Volume analyzed (in voxels)
% V             Image file handles (see spm_vol)
% df            Residual degrees of freedom of raw t-statistic
% MaxT          2xnPerm matrix of [max;min] t-statistics per perm
% ST_Ut         Threshold above which suprathreshold info was collected.
%               Voxel locations, t and perm are saved in SnPM_ST.mat for
%               t's greater than ST_Ut. ST_Ut=Inf if not saving STCdata
%
% Pointlist of (Pseudo) t-statisic for actual labelling:SnPMt.mat
%-----------------------------------------------------------------------
% - saved by spm_snpm
% XYZ		- 1xS matrix of voxel statistics (t or pseudo-t)
%
% Co-ordinates of pointlist:				XYZ.mat
%-----------------------------------------------------------------------
% - saved by spm_snpm
% XYZ		- 3xS matrix of co-ordinates [x;y;z] of voxels on SnPMt
%
% Suprathreshold cluster statistics:			SnPM_ST.mat
%-----------------------------------------------------------------------
% - saved by spm_snpm
% SnPM_ST	- Suprathreshold cluster statistics, see spm_snpm.m
% NB: This file is only required for suprathreshold cluster size analysis



%-Setup
%=======================================================================
fprintf('\nSnPM: spm_snpm_pp\n'),fprintf('%c','='*ones(1,72)),fprintf('\n')

%-Initialise variables & constants
%-----------------------------------------------------------------------
tol = 1e-4;	% Tolerance for comparing real numbers
		% Two reals with abs(a-b)<tol are considered equal
		% ( Reals have to be compared for equality when        )
		% ( computing adjusted p-values                        )

%-SetUp figure window
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph), Fgraph=spm_figure('Create','Graphics'); end
spm_clf(Finter), spm_clf(Fgraph)
set(Finter,'Name','SnPM PostProcess');


%-Get Data
%=======================================================================
% Get analysis directory
if nargin==0
    tmp = spm_get(1,'SnPM.mat','Select SnPM.mat for analysis...');
    CWD  = spm_str_manip(tmp,'hd');
end

%-Load Config file & SnPM permutation data
load(fullfile(CWD,'SnPMcfg'))
load(fullfile(CWD,'SnPM'))


%-Ask whether positive or negative effects be analysed
%-----------------------------------------------------------------------
bNeg = BAT.bNeg; % spm_input('Positive or negative effects?',1,'b','+ve|-ve',[0,1],1);

%-Form full Tmax distribution
%-----------------------------------------------------------------------
%-Tmin are in second column of MaxT, stored with *+ve* values
if bhPerms
	MaxT   = [ MaxT; flipud(fliplr(MaxT)) ];
	PiCond = [PiCond; -flipud(PiCond)];
end
%-Take MaxT for increases or decreases according to bNeg
MaxT = MaxT(:,bNeg+1);
nPerm = size(MaxT,1);
[StMaxT, iStMaxT] = sort(MaxT);

%-Load statistic image
%-----------------------------------------------------------------------
load(fullfile(CWD,'SnPMt'))
load(fullfile(CWD,'XYZ'))

%-Negate if looking at negative contrast
%-----------------------------------------------------------------------
if bNeg
	SnPMt    = -SnPMt;
	CONT     = -CONT;
end

%-Get ORIGIN, etc
DIM    = [V(1).dim(1) V(1).dim(2) V(1).dim(3)];
VOX    = [V(1).mat(1,1) V(1).mat(2,2) V(1).mat(3,3)];
MAT    = V(1).mat;
IMAT   = inv(MAT);
ORIGIN = IMAT(1:3,4);

% Template vol structure
Vs0 = struct('fname',	'',...
	     'dim',	[DIM,spm_type('float')],...
	     'mat',	MAT,...
	     'pinfo',	[1 0 0]',...
	     'descrip',	'');

%-Write out images?
%=======================================================================
%-Write out statistic image
if BAT.writeStat %spm_input('Write out statistic img?','+1','y/n',[1,0],2)
	Fname = 'SPMt';

	%-Don't ask about t2z conversion
	%---------------------------------------------------------------
	bt2z = 0;
	if ~bVarSm
% 	    bt2z = spm_input('Convert t -> z prior to writing?',...
% 	    	'+0','y/n',[1,0],1);
	    tmp = sprintf('SPMt - %d df',df);
	else
	    tmp = 'SnPMt - pseudo t';
	end

	%-Reconstruct statistic image from XYZ & SnPMt
	%---------------------------------------------------------------
	fprintf('Working on statistic image...');
	t = zeros(1,prod(DIM));
	if ~bt2z
		t(spm_xyz2e(XYZ,V)) = SnPMt;
	else
		t(spm_xyz2e(XYZ,V)) = spm_t2z(SnPMt,df);
		tmp = [tmp,' (Gaussianised)'];
	end
	fprintf('done\n')

	%-Write out to analyze file
	%---------------------------------------------------------------
	fprintf('Writing statistic image...');
	if bNeg, Fname = [Fname,'_neg']; tmp = [tmp,' (-ve contrast)']; end
	Vs = Vs0; 
	Vs.fname = Fname; Vs.descrip = tmp;
	Vs = spm_create_image(Vs);
	t = reshape(t,DIM);
	for p=1:Vs.dim(3)
	  Vs = spm_write_plane(Vs,t(:,:,p),p);
	end
	clear t bt2z
	fprintf('done\n')
end

%-Write out filtered statistic image?  (Get's done later)
%-----------------------------------------------------------------------
WrtFlt = 1; %spm_input('Write filtered statistic img?','+1','y/n',[1,0],2);
if WrtFlt
	WrtFltFn = 'SnPMt_filtered';
    if bNeg, WrtFltFn = 'SnPMt_neg_filtered';,end
	%WrtFltFn=spm_input('Filename ?','+1','s',WrtFltFn);
end

%-Write out (full) Single-Step adjusted p-value image?
%-----------------------------------------------------------------------
if 1 % spm_input('Write full SS adj p-value img?','+1','y/n',[1,0],2)

	%-Compute full Single-Step adjusted p-values (watch bNeg)
	%---------------------------------------------------------------
	fprintf('Working on Single Step adjusted p-values...');
	SnPMp = zeros(size(SnPMt));
	for t = MaxT'
		%-Adjusted p is proportion of randomisation greater or
		% equal to statistic.
		%-Use a > b -tol rather than a >= b to avoid comparing
		% two reals for equality.
		SnPMp = SnPMp + (t > SnPMt -tol);
	end
	SnPMp = SnPMp / nPerm;
	fprintf('done\n')

	%-Reconstruct full image from XYZ & p_SSadj
	%---------------------------------------------------------------
	fprintf('Reconstructing Single Step adjusted p-value image...');
	p = zeros(1,prod(DIM));
	p(spm_xyz2e(XYZ,V)) = 1 -SnPMp;
	fprintf('done\n')

	%-Write out to analyze file
	%---------------------------------------------------------------
	fprintf('Writing adjusted p-value image...');
	Fname = 'SnPMp_SSadj';
	tmp = 'SnPMp SingleStep adjusted p-values: 1-p';
	if bNeg, Fname = [Fname,'_neg']; tmp = [tmp,' (-ve contrast)']; end
	Vs = Vs0; 
	Vs.fname = Fname; Vs.descrip = tmp;
	Vs = spm_create_image(Vs);
	p = reshape(p,DIM);
	for s=1:Vs.dim(3)
	  Vs = spm_write_plane(Vs,p(:,:,s),s);
	end
	clear p SnPMp
	fprintf('done\n')
end


%-Get inference parameters
%=======================================================================

%-Get corrected threshold
%-----------------------------------------------------------------------
alpha = BAT.alpha;  %spm_input('Corrected p value for filtering','+1','e',0.05);

%-Compute critical threshold for level alpha test
%-----------------------------------------------------------------------
if alpha < 1;
	c=ceil((1-alpha)*nPerm);
	C_MaxT=StMaxT(c);
else
	%-Just use voxels with +ve valued SnPMt
	C_MaxT=0;
end

%-Ask whether SupraThreshold cluster size test required
%-----------------------------------------------------------------------
%-If chosen alpha specifies a critical threshold less than the threshold
% ST_Ut used to collect suprathreshold data in spm_snpm, then it makes
% no sense to analyse by spatial extent since the voxels are individually
% significant.
%-bST flags whether spatial extent information was collected.)
bSpatEx = bST & exist(fullfile(CWD,'SnPM_ST.mat'))==2;
if bSpatEx, 
    if (C_MaxT > ST_Ut) | (alpha == 1)
	bSpatEx = BAT.bSpatEx;  %spm_input('Assess spatial extent?','+1','y/n',[1,0],2);
    else
	bSpatEx = 0;	
    end
end

%-Get primary threshold for STC analysis if requested
%-----------------------------------------------------------------------
if bSpatEx
    % Save original ST_Ut
    ST_Ut_0 = ST_Ut;
    %-Threshold must be greater or equal to that (ST_Ut) used to collect
    % suprathreshold data in spm_snpm
    %-If a test level alpha has been set, then it there's no sense in having
    % the threshold greater than C_MaxT, above which voxels are individually 
    % significant
    tmp = 0;
    if bVarSm
	%-If using pseudo-statistics then can't use (uncorrected) 
	% upper tail p-values to specify primary threshold
	if alpha==1	% Not filtering on significance
	    while ~(tmp>=ST_Ut)
		tmp = spm_input(sprintf(...
		    'Primary threshold (>%4.2f)',ST_Ut),'+0');
	    end
	else
	    while ~(tmp>=ST_Ut & tmp<C_MaxT)
		sprintf('Threshold (%4.2f<=Ut<%4.2f)',ST_Ut,C_MaxT);
        tmp = BAT.clExtThresh; 
        if BAT.clExtThresh >= C_MaxT, BAT.clExtThresh = C_MaxT * .95;,end
        if BAT.clExtThresh < ST_Ut, BAT.clExtThresh = ST_Ut;,end
        %spm_input(sprintf(...
		    %'Threshold (%4.2f<=Ut<%4.2f)',ST_Ut,C_MaxT),'+0');
	    end
	end
else
	%-Statistic image is t with df degrees of freedom
	pU_ST_Ut  = 1-spm_Tcdf(ST_Ut,df);
	if alpha==1	% Not filtering on significance
	    while ~( tmp>=ST_Ut | (tmp>0 & tmp<=pU_ST_Ut))
		tmp = spm_input(sprintf(...
		    'Threshold (p<=%4.2f I t>=%4.2f)',pU_ST_Ut,ST_Ut),'+0'); 
	    end
	else
	    pU_C_MaxT = 1-spm_Tcdf(C_MaxT,df);
	    while ~((tmp>=ST_Ut & tmp<C_MaxT) | (tmp>pU_C_MaxT & tmp<=pU_ST_Ut))
		tmp = spm_input(sprintf(...
		    'Ut (%4.2f<p<%4.2f I %4.2f<t<%4.2f)',...
		    pU_C_MaxT,pU_ST_Ut,ST_Ut,C_MaxT),'+0');
	    end
	    clear pU_C_MaxT
	end
	clear pU_ST_Ut
	if (tmp < 1), tmp = spm_invTcdf(1-tmp,df); end
    end
    ST_Ut = tmp;
end


%-Show permutation distributions?
%-----------------------------------------------------------------------
%ShwDst = spm_input('Display permutation distribution[s]?','+1','y/n',[1,0],1);
ShwDst = 1;

%=======================================================================
%- C O M P U T A T I O N
%=======================================================================
set(Finter,'Pointer','Watch')

%-Calculate distribution of Maximum Suprathreshold Cluster size
%-Calculate critical Suprathreshold Cluster Size
%=======================================================================
if bSpatEx
	fprintf('Working on spatial extent...\n');

	%-Compute suprathreshold voxels - check there are some
	%---------------------------------------------------------------
	fprintf('\tComputing suprathreshold voxels...');
	Q     = find(SnPMt > ST_Ut);
	SnPMt = SnPMt(Q);
	XYZ   = XYZ(:,Q);
	if isempty(Q)
		set(Finter,'Pointer','Arrow')
		figure(Fgraph)
		axis off
		text(0,0.97,CWD,'Fontsize',16,'FontWeight','Bold');
		text(0,0.93,sprintf('No voxels above threshold %4.2f',ST_Ut));
		ShowDist(MaxT,C_MaxT);
		return
	end
	fprintf('done\n')

	%-Load & condition statistics
	%---------------------------------------------------------------
	fprintf('\tLoading & conditioning SupraThreshold statistics...');
	load(fullfile(CWD,'SnPM_ST'))
	%-SnPM_ST stores columns of [x;y;z;abs(t);perm] with perm negative
	% where the exceedence was t < -ST_Ut_0
	%-Trim statistics according to threshold ST_Ut, if ST_Ut > ST_Ut_0
        if abs(ST_Ut - ST_Ut_0) > 0.01
		tmp = find(SnPM_ST(4,:)>ST_Ut);
		SnPM_ST = SnPM_ST(:,tmp);
		clear tmp;	    
	else
		% We're close enough; we'll just use ST_Ut_0 instead
		ST_Ut = ST_Ut_0;
        end
	%-Negate perm numbers if looking at negative contrast
	if bNeg
		SnPM_ST(5,:) = -SnPM_ST(5,:);
	end
	if bhPerms
		%-Renumber negative perms according to -flipud PiCond
		tQ = SnPM_ST(5,:)<0;
		SnPM_ST(5,tQ) = nPerm +1 +SnPM_ST(5,tQ);
	else
		%-Not bhPerms: Lose entries for negative excursions
		SnPM_ST = SnPM_ST(:,SnPM_ST(5,:)>0);
	end
	fprintf('done\n')

	%-Calculate distribution of Maximum SupraThreshold Cluster size
	%---------------------------------------------------------------
	fprintf('\tComputing dist. of max SupraThreshold cluster size: ');
	MaxSTCS = zeros(nPerm,1);
	SetLvl  = zeros(nPerm,1);
	fprintf('Perms left:     ');
	for i = nPerm:-1:1
                if (rem(i,10)==0); fprintf('\b\b\b\b%-4u',i); end
		tQ = (SnPM_ST(5,:)==i);
		if any(tQ)
         		%-Compute cluster labellings for this perm
			%===== SnPM99 change =================
         		Locs_mm = SnPM_ST(1:3,tQ);
			Locs_mm (4,:) = 1;
			Locs_vox = IMAT * Locs_mm;
			%===== SnPM99 change =================
			tmp = spm_clusters(Locs_vox(1:3,:));
			%-Work out maximum cluster size (honest!)
			MaxSTCS(i) = max(diff(find([diff([0,sort(tmp)]),1])));
			SetLvl(i)  = max(tmp);
		end
	end
	fprintf('\b\b\b\bdone\n');
	%-Save perm 1 stats for use later - [X;Y;Z;T;perm;STCno]
	STCstats = [ SnPM_ST(:,tQ); tmp];

	%-Compute critical SupraThreshold Cluster size
	[StMaxSTCS, iStMaxSTCS] = sort(MaxSTCS);
	if alpha < 1
		C_STCS = StMaxSTCS(c);
	else
		C_STCS = 0;
	end

	%-Check XYZ for points > ST_Ut in perm 1 matches
	% XYZ computed above for SnPMt > ST_Ut
	if ~all(all( SnPM_ST(1:3,SnPM_ST(5,:)==1) == XYZ ))
		error('ST XYZ don''t match between STCS & thresh')
	end
end

%-Save some time consuming results
%-----------------------------------------------------------------------
if bSpatEx, save SnPM_pp STCstats MaxSTCS SetLvl, end


%-Filter data at specified corrected p-value alpha
%=======================================================================
if bSpatEx
    %-Analysing spatial extent
    %-NB:alpha==1 implies C_MaxT==C_STCS==0.
    % Since ST_Ut>0 filtering has no effect if alpha==1, so skip it.
    if alpha<1
	%-Filter on significance of cluster size
	%---------------------------------------------------------------
	fprintf('Filtering on cor.sig. at suprathreshold cluster level...');
	nSTC     = max(STCstats(6,:));
	STCS     = diff(find([diff([0,sort(STCstats(6,:))]),1]));
	Q        = [];
	for i = 1:nSTC
		tQ = find(STCstats(6,:)==i);
		if ( STCS(i) > C_STCS | max(STCstats(4,tQ)) > C_MaxT )
			Q        = [Q tQ];
		end
	end
	if ~isempty(Q)
		SnPMt    = SnPMt(Q);
		XYZ      = XYZ(:,Q);
		STCstats = STCstats(:,Q);
	end
	fprintf('done\n')
    end
else
	%-Truncate at critical threshold for level alpha test
	% NB if alpha==1 then C_MaxT is set to 0, and filter on +ve SnPMt
	fprintf('Filtering on cor.sig. at voxel level...');
	Q     = find(SnPMt > C_MaxT);
	if length(Q)
		SnPMt = SnPMt(Q);
		XYZ   = XYZ(:,Q);
	end
	fprintf('done\n')
end



%-Return if there are no voxels
%-----------------------------------------------------------------------
if isempty(Q)
	set(Finter,'Pointer','Arrow')
	figure(Fgraph)
	axis off
	text(0,0.97,CWD,'Fontsize',16,'FontWeight','Bold');
	tmp='voxels'; if bSpatEx, tmp='suprathreshold clusters'; end
	text(0,0.93,sprintf(...
		'No %s significant at alpha=%6.4f (corrected)',tmp,alpha));
	if bSpatEx,
	  ShowDist(MaxT,C_MaxT,MaxSTCS,C_STCS);
	else	   
	  ShowDist(MaxT,C_MaxT);
	end
	return
end

%-Characterize local excursions in terms of maxima:
% #voxels STC_N; MaxTs STC_SnPMt; locations STC_XYZ, & region# STC_r
%-----------------------------------------------------------------------
%===== SnPM99 change =============================================
TempXYZmm = XYZ;
TempXYZmm(4,:) = 1;
TempXYZvoxel = IMAT*TempXYZmm;
TempXYZvoxel= TempXYZvoxel(1:3,:);

[STC_N, STC_SnPMt, STC_XYZ, STC_r] = spm_max(SnPMt,TempXYZvoxel);

TempXYZvoxel = STC_XYZ;
TempXYZvoxel(4,:) = 1;
TempXYZmm = MAT * TempXYZvoxel;
STC_XYZ = TempXYZmm(1:3,:);
%===== SnPM99 change =============================================

%-Compute adjusted significances for local maxima, & regions (if required)
%-----------------------------------------------------------------------
Pt = ones(size(STC_r));
for i = 1:length(STC_r)
	%-Use a > b -tol rather than a >= b to avoid comparing reals
	Pt(i) = sum(MaxT > STC_SnPMt(i) -tol) / nPerm;
end
if ~bVarSm
	Pu    = 1 - spm_Tcdf(STC_SnPMt,df);
end
if bSpatEx
	%-Compute single step adjusted p-values for region size: pSTSC_SS
	Pn    = ones(size(STC_r));
	for i = 1:length(STC_r)
		Pn(i) = sum(MaxSTCS>=STC_N(i)) / nPerm;
	end
end



%=======================================================================
%-D I S P L A Y
%=======================================================================
figure(Fgraph)


if (ShwDst)
  axis off
  if (bSpatEx)
    text(0,0.97,'Permutation Distribution','Fontsize',16,'FontWeight','Bold');
    ShowDist(MaxT,C_MaxT,MaxSTCS,C_STCS);
  else	     
    text(0,0.97,'Permutation Distributions','Fontsize',16,'FontWeight','Bold');
    ShowDist(MaxT,C_MaxT);
  end
end
spm_print
if BAT.pause, disp('Press <RETURN> to continue'); pause, end

spm_clf(Fgraph)


%-Maximium intenisty projection of SPM{Z}
%=======================================================================
hmip = axes('Position',[0.05 0.5 0.5 0.5]);
spm_mip(SnPMt,XYZ,MAT,DIM); axis image
if bVarSm
	title('SnPM{Pseudo-t}','FontSize',16,'Fontweight','Bold')
else
	title('SnPM{t}','FontSize',16,'Fontweight','Bold')
end

%-Design matrix and contrast
%=======================================================================
hDesMtx = axes('Position',[0.65 0.6 0.2 0.2]);
imagesc((spm_DesMtx('Sca', [H,C,B,G],HCBGnames) + 1)*32)
xlabel 'Design Matrix'
set(hDesMtx,'XTick',[],'XTickLabel','')
hConAxes = axes('Position',[0.65 0.8 0.2 0.1]);
h = bar(CONT(1,:));  hold on
set(h,'FaceColor',[1 1 1]*.8)
tX = get(h,'Xdata'); tY = get(h,'Ydata');
set(gca,'Xlim',[min(tX(:)) max(tX(:))]) 
title 'contrast'; axis off; hold off
 

%-Table of regional effects
%=======================================================================
%-Table headings
%-----------------------------------------------------------------------
hTable = axes('Position',[0.1 0.1 0.8 0.46],...
	'YLim',[0,27],'YLimMode','manual','Visible','off');
y = 26;
text(0,y,['P values & statistics:   ',spm_str_manip(CWD,'a40')],...
	'FontSize',12,'FontWeight','Bold');
y  = y -1;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y = y -0.8;

text(0.00,y,'region','FontSize',10);
text(0.10,y,'size {k}','FontSize',10);
if bSpatEx
	text(0.24,y,'P(K     >= k)','FontSize',10);
	text(0.275,y-0.2,'max','Fontsize',6);
end
if ~bVarSm
	text(0.42,y,'t','FontSize',10);
else
	text(0.42,y,'Pseudo-t','FontSize',10);
end
text(0.5,y,'P( T    >= u)','FontSize',10);
text(0.535,y-0.2,'max','Fontsize',6);
if ~bVarSm, text(0.65,y,'(uncorrected)','Fontsize',8); end
text(0.84,y,'{x,y,z} mm','FontSize',10);
y = y -0.8;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y  = y -1;

%-List of maxima
%-----------------------------------------------------------------------
r = 1;
bUsed = zeros(size(STC_SnPMt));
while max(STC_SnPMt.*(~bUsed)) & (y > 3)

	[null, i] = max(STC_SnPMt.*(~bUsed));	% Largest t value
	j         = find(STC_r == STC_r(i));	% Maxima in same region


	%-Print region and largest maximum
	%-------------------------------------------------------------------
	text(0.00,y,sprintf('%-0.0f',r),...
					'Fontsize',10,'FontWeight','Bold')
	text(0.10,y,sprintf('%-0.0f',STC_N(i)),...
					'Fontsize',10,'FontWeight','Bold')
	if bSpatEx
		text(0.24,y,sprintf('%-0.3f',Pn(i)),...
					'Fontsize',10,'FontWeight','Bold')
	end
	text(0.42,y,sprintf('%-0.2f',STC_SnPMt(i)),...
					'Fontsize',10,'FontWeight','Bold')
	text(0.54,y,sprintf('%-0.3f',Pt(i)),...
					'Fontsize',10,'FontWeight','Bold')
	text(0.82,y,sprintf('%-6.0f',STC_XYZ(:,i)),...
					'Fontsize',10,'FontWeight','Bold')
	if ~bVarSm
		text(0.65,y,sprintf('(%-0.6f)',Pu(i)),'Fontsize',10)
	end
	y = y -1;

	%-Print up to 3 secondary maxima (>8mm apart)
	%-------------------------------------------------------------------
	[null, k] = sort(-STC_SnPMt(j));	% Sort on t value
	D         = i;
	for i = 1:length(k)
	    d     = j(k(i));
	    if min( sqrt( sum((STC_XYZ(:,D) - ...
			STC_XYZ(:,d)*ones(1,size(D,2))).^2) ) ) > 8;
		if length(D) < 3
		    text(0.42,y,sprintf('%-0.2f',STC_SnPMt(d)),'Fontsize',10)
		    text(0.54,y,sprintf('%-0.3f',Pt(d))       ,'Fontsize',10)
		    text(0.82,y,sprintf('%-6.0f',STC_XYZ(:,d)),'Fontsize',10)
		    if ~bVarSm
			text(0.65,y,sprintf('(%-0.6f)',Pu(d)), 'Fontsize',10)
		    end
		    D = [D d];
		    y = y -1;
		end
	    end
	end

	bUsed(j) = (bUsed(j) | 1 );		%-Mark maxima as "used"
	r = r + 1;				% Next region
end
clear i j k D d r


%-Footnote with SnPM parameters
%=======================================================================
line([0,1],[0.5,0.5],'LineWidth',1,'Color',[0 0 0])
y = 0;
if bSpatEx
    tmp = sprintf('Threshold = %7.4f',ST_Ut);
    if ~bVarSm
	tmp=[tmp,sprintf(' (p = %6.4f)',spm_Tcdf(-ST_Ut,df))];
    end
    text(0,y,tmp,'FontSize',8)
    text(0.7,y,sprintf('Critical STCS = %d voxels',C_STCS),'FontSize',8)
    y = y -0.8;
end
text(0,y,sprintf('alpha = %6.4f, df = %d',alpha,df),'FontSize',8)
text(0.7,y,sprintf('Critical threshold = %7.4f',C_MaxT),'FontSize',8)
y = y -0.8;
text(0,y,sprintf('Volume = %d %5.2fx%5.2fx%5.2f mm voxels',...
	S,VOX(1),VOX(2),VOX(3)),'FontSize',8);
y = y -0.8;
text(0,y,sprintf('Design: %s',sDesign),'FontSize',8);
y = y -0.8;
text(0,y,sprintf('Perms: %s',sPiCond),'FontSize',8);
if bVarSm
	y = y -0.8;
	text(0,y,sVarSm,'FontSize',8)
end

spm_print

set(Finter,'Pointer','Arrow')

%- Image output?
%=======================================================================
%-Write out filtered SnPMt?
if WrtFlt
  
	Fname = WrtFltFn;

	%-Dont ask about t2z conversion
	%---------------------------------------------------------------
	bt2z = 0;
	if ~bVarSm
% 	    bt2z = spm_input('Convert t -> z prior to writing?',...
% 	    	'+1','y/n')=='y';
	    tmp = sprintf('SPMt - %d df',df);
	else
	    tmp = 'SnPMt - pseudo t';
	end

	%-Reconstruct filtered image from XYZ & SnPMt
	%---------------------------------------------------------------
	t = zeros(1,prod(DIM));
	if ~bt2z
		t(spm_xyz2e(XYZ,V)) = SnPMt;
	else
		t(spm_xyz2e(XYZ,V)) = spm_t2z(SnPMt,df);
		tmp = [tmp,' (Gaussianised)'];
	end
	if ~bSpatEx
		tmp=sprintf('%s p<%10g corrected @ voxel level',tmp,alpha)
	elseif bt2z
		tmp=sprintf('%s p<%10g corrected @ cluster level, u=%4.2',...
	    		tmp,alpha,spm_t2z(ST_Ut,df));
	else
		tmp=sprintf('%s p<%10g corrected @ cluster level, u=%4.2',...
	    		tmp,alpha,ST_Ut,df);
	end

	%-Write out to analyze file
	%---------------------------------------------------------------
	Vs = Vs0; 
	Vs.fname = Fname; Vs.descrip = tmp;
	Vs = spm_create_image(Vs);
	t = reshape(t,DIM);
	for p=1:Vs.dim(3)
	  Vs = spm_write_plane(Vs,t(:,:,p),p);
	end
	clear t
end

%-Reset Interactive Window
%-----------------------------------------------------------------------
spm_figure('Clear','Interactive')



function ShowDist(T,cT,C,cC)
%
% Display permutation distributions on current figure, using position
% pos.
%
% We assume that the observed (aka correctly labeled data) are the first
% element in the distribution.
%

pos1 = [0.125 0.50 0.75 0.3];
pos2 = [0.125 0.08 0.75 0.3];

%
% Display intensity perm dist
%
% Bin width rule from Scott, "Multivariate Density Estimation", 1992, pg 55.
%
axes('position',pos1)
BinWd  = 3.5*std(T)*length(T)^(-1/3);  
nBin   = floor((max(T)-min(T))/BinWd)+1;
nBin   = min(max(nBin,10),50);

hist(T,nBin);
set(get(gca,'Children'),'FaceColor',[.5 .5 .5]) 
title('Permutation Distribution:  Maximum Statistic','FontSize',14)
Ylim = get(gca,'ylim'); Xlim = get(gca,'xlim');
line(T(1)*[1 1],Ylim.*[1 0.95],'LineStyle',':');
text(T(1)+diff(Xlim)*0.01,Ylim(2)*0.95,'Observed','FontSize',10)
line(cT*[1 1],Ylim.*[1 0.85],'LineStyle','-');
text(cT+diff(Xlim)*0.01,Ylim(2)*0.85,'Threshold','FontSize',10);

if (nargin>2)

  %
  % Display cluster size perm dist
  %

  axes('position',pos2)
  BinWd  = 3.5*std(C)*length(C)^(-1/3);
  nBin   = floor((max(C)-min(C))/BinWd)+1;
  nBin   = min(max(nBin,10),50);

  hist(C,nBin);
  set(get(gca,'Children'),'FaceColor',[.5 .5 .5]) 
  title('Permutation Distribution: Maximum Cluster Size','FontSize',14)
  Ylim = get(gca,'ylim'); Xlim = get(gca,'xlim');
  line(C(1)*[1 1],Ylim.*[1 0.95],'LineStyle',':');
  text(C(1)+diff(Xlim)*0.01,Ylim(2)*0.95,'Observed','FontSize',10)
  line(cC*[1 1],Ylim.*[1 0.85],'LineStyle','-');
  text(cC+diff(Xlim)*0.01,Ylim(2)*0.85,'Threshold','FontSize',10)

end

return
