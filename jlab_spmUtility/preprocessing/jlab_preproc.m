function jlab_preproc(OPT)

spm fmri 	%this should start spm99


% check to make sure subject directory exists
% --------------------------------------------------------------------
str = [OPT.studydir filesep OPT.subjcode];
if ~(exist(str) == 7), error([str ' does not exist!']), end

% check for required parameters
% --------------------------------------------------------------------
req = {'studydir' 'subjcode' 'nruns' 'coreg_anat' 'target' ...
	'segment_object' 'unseg_object' 'object' 'canonicalT1' ...
	'createmean' 'detNorm' 'applyNorm' 'smooth' 'voxsize' ...
	'kernel', 'checknorm', 'bigmask'};
for i = 1:length(req)
	if ~isfield(OPT,req{i}), 
		error(['jlab_preproc -> missing required field:' req{i}])
	end
end

% change to directory and set up some preliminary variables
% --------------------------------------------------------------------
eval(['cd ' str])

subjDir = str;
t1dir = [str filesep 'anatomy'];


% --------------------------------------------------------------------
% * make mean functional for each subject
% --------------------------------------------------------------------
if OPT.createmean

   % create mean image in scan1 directory
   % all images should be in separate dirs for each run, labeled scan1, scan2,...etc.

   disp('Computing mean from all imgs in scan* dirs')
   for i = 1:OPT.nruns
   	eval(['cd ' subjDir])
   	myP = getfiles(['scan' num2str(i) filesep 'ra*.img']);
  	disp(['	scan' num2str(i) ': found ' num2str(length(myP)) ' images....computing...'])
   	myP = str2mat(myP);
   	tor_spm_mean_ui(myP)

   	!rm mean.mat
   	eval(['!mv mean.img mean_scan' num2str(i) '.img'])
   	eval(['!mv mean.hdr mean_scan' num2str(i) '.hdr'])
   end
   eval(['cd ' subjDir])
end

% --------------------------------------------------------------------
% * coregister
% --------------------------------------------------------------------	

if OPT.coreg_anat

	eval(['cd ' t1dir]) 	% get to anatomy directory
	
	%spm_coregister(target image,object image,affnorm target to,affnorm object to, others, flags)
	%middle 2 are same modality templates for coregistration btwn modalities e.g. T1 to T2

	if OPT.coreg_2_funct
		spm_coregister(OPT.target, OPT.object)
	else
		spm_coregister(OPT.target, OPT.object) 	% you MUST be in SPM99 "fMRI time-series" to do this	
										% otherwise it can't display: "Error using ==> axes"
	end
	pause(4)
	spm_print

	% 1/26/02 do not reslice
	% spm_reslice(myt1,struct('mask',0,'mean',0,'hold',1,'which',2));
    eval(['cd ' subjDir])
end

if OPT.check_coreg
    mytarget = OPT.target;	% for functionals w full path name.
    %[d f e] = fileparts(OPT.target); 
    %mytarget = fullfile('anatomy',[f e]);  % ('anatomy',[f e]);  %('scan1',[f e]);
    myobject = fullfile('anatomy',deblank(OPT.object));
    P = str2mat(mytarget,myobject);
  	spm_check_registration(P);
   	spm_orthviews('Interp',0);
	spm_orthviews('window',1,[0 6000]); spm_orthviews('window',2,[0 250])
   	title(['Coreg ' OPT.object ' to ' OPT.target ' subject ' OPT.subjcode])
    input('Surf through image and press return when ready to print.')
        
    spm_print			%this prints the results of the coregistration to spm.ps file
    
    % for reorienting by hand...
    if OPT.handadjust
    	disp('Make sure the images line up by surfing through them.')
    	disp('Adjust by using the trans/rotation factors in spm_image')
    	disp(['spm_image(''init'',myobject)'])
    	disp(['to check again, use:'])
    	disp('spm_check_registration(P)')
    	disp(['Your images are:'])
    	P
        disp(['To adjust intensity: spm_orthviews(''window'',1,[0 1500])'])
	disp(['When done aligning: spm_coregister(mytarget, myobject)'])
    	spm_image('init',myobject)
    	keyboard
   end	
end


% --------------------------------------------------------------------
% * segment / print for each subject
% --------------------------------------------------------------------

if OPT.segment_object

	% FORMAT spm_segment(PF,PG,opts)
	% PF   - name(s) of image(s) to segment (must have same dimensions).
	% PG   - name(s) of template image(s) for realignment.
	%      - or a 4x4 transformation matrix which maps from the image to
	%        the set of templates.
	% opts - options string.
	%        - 't' - write images called *_seg_tmp* rather than *_seg*
	%                (that are smoothed with an 8mm Gaussian).
	%        - 'f' - fix number of voxels in each cluster
	%        - 'c' - attempt to correct small intensity inhomogeneities
	%        - 'C' - attempt to correct large intensity inhomogeneities
	%        - 'w' - write inhomogeneity corrected image(s) (with the 'c'
	%                or 'C' options only)

	spm_segment(fullfile(t1dir,OPT.unseg_object),OPT.canonicalT1,'C')
	% spm_print
	
end


if OPT.detNorm | OPT.applyNorm | OPT.smooth | OPT.twostagenorm
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	% Obligatory first part of the automation script                           % 
	%*************************************************************************** 
	%  fmriTEST = 1 to run the simulation mode 
	%           = 0 to run the actual post-processing 
	%  fmriDIR  = the directory in which all your anatomical data is situated 
	%--------------------------------------------------------------------------- 
	global fmriTEST fmriDIR;    %set global variables
	fmriDIR = OPT.studydir; 	%this is the directory where all the expt data is
	fmriTEST = 0; 
	f99_CWD('LOG');             %make the LOG and RESULTS directory 
	eval(['cd ' fmriDIR])

	% must be running in spm, otherwise 
	% starting orientations - pick one in function call later.
	% if you use -R in reconstruction, you're in Neuro - as of 8/23/01
	neuro = [0 0 0 0 0 0 1 1 1 0 0 0];
	radio = [0 0 0 0 0 0 -1 1 1 0 0 0];

	if ~isfield(OPT,'twostagenorm'), OPT.twostagenorm = 0;, end

	objectloc = [OPT.subjcode filesep 'anatomy' filesep OPT.object];
		%set objectloc so that we get to the correct directory where the object T1 is

	
	if OPT.twostagenorm
		secondobjectloc = [OPT.subjcode filesep 'anatomy' filesep 'n' OPT.object];
	end

	if OPT.applyNorm | OPT.smooth
		% check for number of scans found - and adjust if necessary!
		% ---------------------------------------------------------------------
		eval(['D=dir(''' OPT.subjcode ''');'])
		m = strmatch('scan',str2mat(D.name));
		if length(m) ~= OPT.nruns
			warning(['Only ' num2str(length(m)) ' scan directories found.  Adjusting OPT.nruns.'])
			OPT.nruns = length(m);
		end
	end

end

if OPT.detNorm

	% --------------------------------------------------------------------
	% * determine normalization parameters
	% --------------------------------------------------------------------

 	f99_sn3d(1,...							% option 1 = 'Determine Parameters Only' 
 	objectloc, ...							% determine params from
 	0, ...									% masking
	'',... 									% name of mask image
 	[],... 	% directories for locating images...relative path to fmriDIR
 	[],... 									% EXP: file name wildcards for imgs to normalize
 	OPT.canonicalT1, ...					% template
 	0, ...									% mask template?
	[],... 									% name of mask img for template			
 	'affine3', ...							% type of affine reg?
 	neuro, ...								% neuro or radio convention
 	13, ...									% number of basis functions			
 	16, ...									% nonlinear iterations
	0.01,... 								% nonlinear regularization [.01 default]					
        1, ...									% interpolation method (1 = bilinear, default)
	[-78 78 -112 76 -50 85], ...			% bounding box
	OPT.voxsize, ...						% voxel size
	'');									% other params & voxel size - bilinear interp
	
	spm_print
end

if OPT.applyNorm

	% --------------------------------------------------------------------
	% * build list of where all ra_vols are and dir wildcards
	% --------------------------------------------------------------------
	% this will crash if you have more than 9 scans.	

	for i = 1:OPT.nruns		
		ra_dir(i,:) = [OPT.subjcode filesep 'scan' num2str(i)];
		exp(i,:) = 'ravol*img';
	end

	% --------------------------------------------------------------------
	% * apply normalization to all ra_vols
	% --------------------------------------------------------------------

 	f99_sn3d(2,...									% option 2 = 'Write normalized Only' 
 		objectloc, ...							% determine params from
 		0, ...									% masking
		'',... 									% name of mask image
 		[ra_dir],... 							% directories for locating images
 		[exp],... 								% EXP: file name wildcards for imgs to normalize
 		[], ...									% template
 		0, ...									% mask template?
		[],... 									% name of mask img for template			
 		'affine3', ...							% type of affine reg?
 		neuro, ...								% neuro or radio convention
 		13, ...									% number of basis functions			
 		16, ...									% nonlinear iterations
		0.01,... 								% nonlinear regularization [.01 default]					
       		1, ...									% interpolation method (1 = bilinear, default)
		[-78 78 -112 76 -50 85], ...			% bounding box
		OPT.voxsize, ...						% voxel size
		'');									% other params & voxel size - bilinear interp

end

if OPT.twostagenorm

	% --------------------------------------------------------------------
	% * determine normalization parameters
	% --------------------------------------------------------------------

 	f99_sn3d(1,...							% option 1 = 'Determine Parameters Only' 
 	secondobjectloc, ...					% determine params from
 	0, ...									% masking
	'',... 									% name of mask image
 	[],... 	% directories for locating images...relative path to fmriDIR
 	[],... 									% EXP: file name wildcards for imgs to normalize
 	OPT.secondcanonicalT1, ...				% template
 	0, ...									% mask template?
	[],... 									% name of mask img for template			
 	'affine3', ...							% type of affine reg?
 	neuro, ...								% neuro or radio convention
 	13, ...									% number of basis functions			
 	16, ...									% nonlinear iterations
	0.01,... 								% nonlinear regularization [.01 default]					
        1, ...									% interpolation method (1 = bilinear, default)
	[-78 78 -112 76 -50 85], ...			% bounding box
	OPT.voxsize, ...						% voxel size
	'');									% other params & voxel size - bilinear interp

	% --------------------------------------------------------------------
	% * build list of where all ra_vols are and dir wildcards
	% --------------------------------------------------------------------
	% this will crash if you have more than 9 scans.	

	clear ra_dir, clear exp
	for i = 1:OPT.nruns		
		ra_dir(i,:) = [OPT.subjcode filesep 'scan' num2str(i)];
		exp(i,:) = 'nravol*img';
	end

	% --------------------------------------------------------------------
	% * apply normalization to all nra_vols
	% --------------------------------------------------------------------

 	f99_sn3d(2,...									% option 2 = 'Write normalized Only' 
 		secondobjectloc, ...					% determine params from
 		0, ...									% masking
		'',... 									% name of mask image
 		[ra_dir],... 							% directories for locating images
 		[exp],... 								% EXP: file name wildcards for imgs to normalize
 		[], ...									% template
 		0, ...									% mask template?
		[],... 									% name of mask img for template			
 		'affine3', ...							% type of affine reg?
 		neuro, ...								% neuro or radio convention
 		13, ...									% number of basis functions			
 		16, ...									% nonlinear iterations
		0.01,... 								% nonlinear regularization [.01 default]					
       		1, ...									% interpolation method (1 = bilinear, default)
		[-78 78 -112 76 -50 85], ...			% bounding box
		OPT.voxsize, ...						% voxel size
		'');									% other params & voxel size - bilinear interp
	
	% remove old nravols to save space! - better to do manually.
	% eval(['!rm ' OPT.subjcode filesep 'scan*' filesep 'nravol*'])

end

	

if OPT.smooth

	if OPT.twostagenorm
		myexp = ['nnravol*.img'];
	else
		myexp = ['nravol*.img'];
	end

	% --------------------------------------------------------------------
	% * smoothing
	% --------------------------------------------------------------------
	for i = 1:OPT.nruns 
  
		f99_smooth(OPT.kernel,...
			[OPT.subjcode filesep 'scan' num2str(i)],... 
			myexp);
	end

end

eval(['cd ' subjDir])

if OPT.checknorm

	% --------------------------------------------------------------------
	% * check normalization parameters
	% --------------------------------------------------------------------

	um_check_norm(OPT)

end

eval(['cd ' subjDir])
%eval(['!mv spm99.ps ' OPT.subjcode '.ps'])



if OPT.bigmask
    
    TOR.in_create = 'y';		% create bigmask image file
    TOR.in_segment = 'n';		% use already created segmented images
    TOR.in_modify = 'n';		% modify spmCFG.mat before estimating
    TOR.in_estim = 'n';			% estimate after modifying spmCFG.mat
    TOR.PG = OPT.canonicalT1;
    TOR.PF = fullfile(subjDir,'anatomy',['n' OPT.object]);	

	% image to segment

    if ~(exist(TOR.PF) == 2), 
        warning([TOR.PF ' does not exist.  Cannot create bigmask.img'])
    else
        tor_glm_specmask(TOR)
    end
end