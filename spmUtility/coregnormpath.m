% 9/11/01
% WARNING:
% 	Make sure origins are set correctly for anatomical and mean functional images!
%	Anatomicals should be set before segmentation.
%	Anatomicals and functionals should be set before coregistering and normalization.
%
% ALSO:
%	Make sure the SPM gui is running before you run this script and normsmooth_modelX
%	Sometimes I get an error in segmenting - needs a string input and doesn't have one.
%	Exiting and restarting spm fmri seems to fix this.

% --------------------------------------------------------------------
% * which operations to perform user options.
% --------------------------------------------------------------------

% anatomical options
	homocorrect = 1;		% Kalina Christoff / Gary Glover's homogeneity correction for spiral images.
	docoregister = 0;
	dosegment = 1;			

% functional options
	createmean = 0;			% create mean functional image - spm's is shifted over for me.

% flow options
	runnormsmooth = 1;		% continue and run all normalization

% substyle					% choose 1 for numbers, 2 for letter/number codes
	substyle = 2;

t1name = 't1gre.img';			% name of original t1 image to use for normalization, etc.

% --------------------------------------------------------------------
% * path and template user options
% --------------------------------------------------------------------

% for jonides7
% -------------------------------------------------------------
%fmriDIR  = '/data4/intext2'; 
%eval(['cd ' fmriDIR])

% normalizeation template(s)
canonicalT1 = '/data3/biman/spm99/spm99/templates/T1.img';
% myavg = '/data1/intext/avganat/avgt1.img';
canonicalfunct = '/data3/biman/spm99/spm99/templates/EPI.img';
segtemplate = '/data3/biman/spm99/spm99/templates/T1_seg1.img';

% for picasso
% -------------------------------------------------------------
fmriDIR  = '/data/biman5'; 
eval(['cd ' fmriDIR])

% normalizeation template(s)
canonicalT1 = '/usr/private/spm99/templates/T1.img';
canonicalfunct = '/usr/private/spm99/templates/EPI.img';
segtemplate = '/usr/private/spm99/templates/T1_seg1.img';



% prompt for subs index if it doesn't exist
%if ~(exist(subs{1}) == 1),
%	substyle = input('Subject list in numbers (1) or letter/number strings (2)? ');
%	if substyle == 1
%		subs = input('Type in string of subject ID numbers: ');
%	elseif substyle == 2
%		subs = input('Type in cell vector of subject ID codes, with {} around it and single quotes: ','s');
%	else
%		error('Uknown substyle - choose 1 or 2.')
%	end
%end


for JJ = subs
	if substyle == 2
   		% if subs get letter/number code   
   		snum = JJ{1};
	elseif substyle == 1
   		% if subs are numbered
   		snum = ['sub' num2str(JJ)];
   	else
		error('Unknown substyle - choose 1 or 2')
	end

   disp(['* Subject ' snum ' ______________________________________________________________'])


   % exit if subject doesn't exist
   str = [fmriDIR filesep snum filesep 'task'];
   if ~(exist(str) == 7), error([str ' does not exist!']), end

   % create task directory if doesn't exist
   str = [fmriDIR filesep snum filesep 'task'];
   if ~(exist(str) == 7), eval(['!mkdir ' str]), end
  
   % --------------------------------------------------------------------
   % * individual subject user options
   % --------------------------------------------------------------------

   % image to determine params from
   % flip all t1's in the x direction, as of 8/23/01, for UM analysis
   % check whether this is needed by comparing anatomicals and functionals.
	% ORIGINAL ANATOMICAL BEFORE SCRIPT PROCESSING

        t1path = [snum '/anatomy'];
	myt1 = [t1path filesep t1name]; 
   

   % coregisters without reslicing, therefore coreg is ft1.img (or your name in myt1)
   % coregistration info is in .mat file of ft1.img (or name in myt1)
   % but then reslices...so it's rft1, or if homogeneity corrected, rhft1.
	% don't use the next line, because myt1 is updated throughout the script:
   	% myrt1 = ['sub' snum '/anatomy/rhft1.img'];

   % --------------------------------------------------------------------
   % * homogeneity correction
   % --------------------------------------------------------------------
if homocorrect
   
   ht1name = ['h' t1name];
   vol_homocor2(myt1,[t1path filesep ht1name]);
   myt1 = [t1path filesep ht1name];
   
   P = str2mat([t1path filesep t1name],myt1);
   spm_check_registration(P);
   spm_orthviews('Interp',0);
   title(['Pre- and post-Homogeneity correction, subject ' snum])
   spm_print

end

   % --------------------------------------------------------------------
   % * make mean functional for each subject
   % --------------------------------------------------------------------
if createmean

   % create mean image in scan1 directory
   % all images should be in separate dirs for each run, labeled scan1, scan2,...etc.

   disp('Computing mean from all imgs in scan* dirs')
   eval(['cd ' fmriDIR filesep snum filesep 'task/'])
   myP = getfiles(['scan1/ra*.img']);
   disp(['	Found ' num2str(length(myP)) ' images.'])
   myP = str2mat(myP);
   tor_spm_mean_ui(myP)

   !rm mean.mat
   disp('Choose mean image and set origin.')
   spm_image
   keyboard
   
   eval(['cd ' fmriDIR])
   
end

   % --------------------------------------------------------------------
   % * get the name of the first and mean images
   % --------------------------------------------------------------------
   disp('Finding the first functional image')
   myP = getfiles([snum '/task/scan1/ra*0001.img'])
   firstfun = myP{1};
   
   meanfun = [snum '/task/mean.img'];


   % --------------------------------------------------------------------
   % * coregister / print / segment / print for each subject
   % --------------------------------------------------------------------	

if docoregister
	%spm_coregister(tgt-template,obj,affnorm target to,affnorm object to, others)
	% last 2 are same modality templates for coregistration btwn modalities.

	spm_coregister(meanfun,myt1,canonicalfunct,canonicalT1,[],[])
	spm_print

	spm_reslice(myt1,struct('mask',0,'mean',0,'hold',1,'which',2));
	rt1name = ['r' ht1name];
	myt1 = [t1path filesep rt1name];
else
	rt1name = [ht1name];
	myt1 = [t1path filesep rt1name];
end

   % --------------------------------------------------------------------
   % * segment / print for each subject
   % --------------------------------------------------------------------

if dosegment

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

	spm_segment(myt1,canonicalT1,'C')
	spm_print

end


end % loop through subjects


% ...now run normsmooth_modelX
% with seg_ template and t1 options.

if runnormsmooth
	normsmooth_modelX
end
