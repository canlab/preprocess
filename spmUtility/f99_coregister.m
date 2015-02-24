function f99_coregister(option,Use_MutInf,respt,targetima,respo,objectima,DIR,EXP)

% '========================================================================='
%    'COREGISTRATION'
% '========================================================================='
% '               --- The Script Parameters Explained ---                   '


% --------------------------------------------------------------------------
% 'do what?'
%    1 - coregister only
%    2 - reslice only
%    3 - coregister & reslice
% 'Use Mutual Information?'
%    1 = yes ; 0 = no
% 'modality of target image'
%    1 = PET | 2 = T1 MRI | 3 = T2 MRI 
%    4 = PD MRI | 5=  EPI | 6 = Transm
% 'target image'
%    e.g. 'run1/first_epi.img'
% 'modality of object image'
%    see target
% 'object image'
%    e.g. 'anat/mpr.img'
% 'Other images to reslice'
% 'DIR'   
%    [matrix containing directories of experiments0] 
% 'EXP'
%    [matrix containing directories of images wildcard] 
% --------------------------------------------------------------------------


%function spm_coregister(opt)
% Between and within mode image coregistration.
%
%____________________________________________________________________________
%
% The TARGET image is the image to which the OBJECT image is realigned.
% If there are any OTHER images, then the same transformations are applied to
% these images as are applied to the OBJECT image.
%
% eg 1) to realign a structural MR image to a sequence of PET images:
%  TARGET: meanPET1.img
%  OBJECT: MRI.img
%   OTHER: -
%
% eg 2) to realign a sequence of PET images to a structural MR image:
%  TARGET: MRI.img
%  OBJECT: meanPET1.img
%   OTHER: PET1.img PET2.img PET3.img etc...
%____________________________________________________________________________
%
% Normally the program has two modes of operation:
% 1) If the modalities of the target image(s) and the object image(s) are
%    the same, then the program performs within mode coregistration by
%    minimising the sum of squares difference between the target and object.
%
% 2) If the modalities differ, then the following is performed:
%    i)   Affine normalisation of object to a template of the same modality,
%         and affine normalisation of the target to a template of the same
%         modality.  Only the parameters which describe rigid body
%         transformations are allowed to differ between these normalisations.
%         This produces a rough coregistration of the images.
%    ii)  The images are partitioned into gray matter, white matter, csf and
%         (possibly) scalp using spm_segment.m.  The mappings from images to
%         templates derived from the previous step are used to map from the
%         images to a set of a-priori probability images of GM, WM and CSF.
%    iii) These partitions are then registered together simultaneously, using
%         the results of step i as a starting estimate.
%
% It is also possible to use Mutual Information image registration.  This
% option is available via the <Defaults> button.
% The Mutual Information is essentially given by:
% H  = H/(sum(H(:))+eps);
% s1 = sum(H,1);
% s2 = sum(H,2);
% H  = H.*log2((H+eps)./(s2*s1+eps));
% mi = sum(H(:));
% where H is a 256x256 histogram, and mi is the mutual information.
% 
% Realignment parameters are stored in the ".mat" files of the "object" and
% the "other" images.
%____________________________________________________________________________
% Refs:
%
% Ashburner J & Friston KJ (1997) Multimodal Image Coregistration and
% Partitioning - a Unified Framework. NeuroImage 6:209-217
%
% Ashburner J, Neelin P, Collins DL, Evans AC & Friston KJ (1997)
% Incorporating Prior Knowledge into Image Registration. NeuroImage 6:344-352
%
% For MI method:
% A Collignon, F Maes, D Delaere, D Vandermeulen, P Suetens & G Marchal
% (1995) "Automated Multi-modality Image Registration Based On
% Information Theory". In the proceedings of Information Processing in
% Medical Imaging (1995).  Y. Bizais et al. (eds.).  Kluwer Academic
% Publishers.
%____________________________________________________________________________
% @(#)spm_coreg_ui.m	2.1 John Ashburner FIL 99/10/13

if nargin~=0,
	if strcmp(lower(option),'defaults'),
		edit_defs;
		return;
	end;
end;

global sptl_UsMtlInfrmtn

global SWD BCH fmriTEST fmriDIR logfile;
SWD = spm('Dir');
DIR1 = [fullfile(SWD,'templates'),filesep];

% fmri init stuff
%----------------
if nargin ~= 8
   error('This script needs 7 arguments');
end

D = [];
for i = 1:size(DIR,1)
   d = [deblank(fmriDIR) filesep deblank(DIR(i,:))];
   D = [D; d];
end
DIR = D;

if f99_exist(fmriDIR,'RESULTS') == 0
   str = ['!mkdir ' fmriDIR  filesep 'RESULTS'];
   eval(str);
end
if f99_exist(fmriDIR,'LOG') == 0
   str = ['!mkdir ' fmriDIR filesep 'LOG'];
   eval(str);
end
if f99_exist([fmriDIR filesep 'LOG'],'COREGISTER.OK')
   str = ['!rm -f ' fmriDIR filesep 'LOG' filesep 'COREGISTER.OK'];
   eval(str);
end

% Start logging
%--------------
logfile = [fmriDIR filesep 'LOG' filesep 'coregister.log'];
tijd=spm('time');
disp('');
disp('*****************************************************************');
f99_log(logfile,['Coregister Script started at ' mat2str(tijd)]);
str = ['Start logging in ' logfile];
f99_log(logfile,str);
%------------------------------------------------------------------------

% mutual information??
% --------------------
%if any(sptl_UsMtlInfrmtn),
if any(sptl_UsMtlInfrmtn) | Use_MutInf
%	spm_mireg_ui;
        f99_mireg(option,respt,targetima,respo,objectima,DIR,EXP);
	return;
end;


%SPMid = spm('FnBanner',mfilename,'2.1');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Coregister');
spm_help('!ContextHelp','spm_coreg_ui.m');

% get number of subjects
%nsubjects = spm_input('number of subjects',1, 'e', 1,...
%	'batch',{},'subject_nb');
nsubjects =1;

if (nsubjects < 1)
	spm_figure('Clear','Interactive');
	return;
end

%p = spm_input('Which option?',2,'m',...
%	'Coregister only|Reslice Only|Coregister & Reslice',...
%	[1 2 3],3, 'batch',{},'opt');
p = option;    
    
if p == 1 | p == 3,
        f99_log(logfile,'Coregistering: ');
	flags = ' ';
	templates = str2mat([DIR1 'PET.img'], ...
		[DIR1 'T1.img'], [DIR1 'T2.img'],...
		[DIR1 'PD.img'], [DIR1 'EPI.img'],...
		[DIR1 'Transm.img'],[DIR1 'SPECT.img']);

	% Get modality of target
	%-----------------------------------------------------------------------
%	respt = spm_input('Modality of first target image?',3,'m',...
%		['target - PET|target - T1 MRI|target - T2 MRI|target - PD MRI|'...
%		 'target - EPI|target - Transm|target - SPECT'],...
%		[1 2 3 4 5 6 7],1,'batch',{},'target_mod');
        strp = ['PET   ';...
                'T1 MRI';...
                'T2 MRI';...
                'PD MRI';...
                'EPI   ';...
                'Transm'];
	PGG = deblank(templates(respt,:));

	% Get modality of object
	%-----------------------------------------------------------------------
%	respo = spm_input('Modality of first object image?',4,'m',...
%		['object - PET|object - T1 MRI|object - T2 MRI|object - PD MRI|'...
%		 'object - EPI|object - Transm|object - SPECT'],...
%		[1 2 3 4 5 6 7],2,'batch',{},'object_mod');
	PFG = deblank(templates(respo,:));

	if respo==6 | respt==6,
		% only perform the first step of the registration
		% because transmission/CT images do not segment very
		% well.
		flags = [flags 'n'];
	end;
	if respt == respo,
		n_images = 1;
	else,
		n_images = Inf;
	end;

	for i = 1:nsubjects,
		% select target(s)
%		if isempty(BCH),
%			PGF = [];
%			while size(PGF,1)<1,
%				PGF = spm_get(n_images,'.img',...
%					['select target image for subject ' num2str(i)]);
%			end;
%		else,
%			PGF = spm_input('batch',{},'target_image');
%		end;
                PGF = [];
                PGF = [fmriDIR filesep targetima];
                f99_log(logfile,[' The target image is ',PGF]);
                f99_log(logfile,['   and this is a ',strp(respt,:)]);

		% select object(s)
%		if isempty(BCH),
%			PFF = [];
%			while size(PFF,1)<1,
%				PFF = spm_get(n_images,'.img',...
%					['select object image for subject ' num2str(i)]);
%			end;
%		else,
%			PFF = spm_input('batch',{},'object_image');
%		end;
                PFF = [];
                PFF = [fmriDIR filesep objectima];
                f99_log(logfile,[' The object image is ',PFF]);
                f99_log(logfile,['   and this is a ',strp(respo,:)]);


		% select others
%		if isempty(BCH),
%			others = spm_get(Inf,'.img',...
%				['select other images for subject ' num2str(i)]);
%		else,
%			PFF = spm_input('batch',{},'other_image');
%		end;
                others = f99_P(DIR,EXP);

		eval(['PGF'    num2str(i) ' = PGF;']);
		eval(['PFF'    num2str(i) ' = PFF;']);
		eval(['others' num2str(i) ' = others;']);
	end;
end;

if p==2,
	for i = 1:nsubjects,
		% select target space
%		if isempty(BCH),
%			PGF = spm_get(1,'.img',...
%				['select image defining space for subject ' num2str(i)]);
%		else,
%			PGF = spm_input('batch',{},'target_image');
%		end;
                PGF = [fmriDIR filesep targetima];

		% select images to reslice
%		if isempty(BCH),
%			PFF = [];
%			PFF = spm_get(Inf,'.img',...
%				['select images to reslice ' num2str(i)]);
%		else,
%			PFF = spm_input('batch',{},'reslice_image');
%		end;
                PFF = f99_P(DIR,EXP);

		eval(['PGF'    num2str(i) ' = PGF;']);
		eval(['PFF'    num2str(i) ' = PFF;']);
		eval(['others' num2str(i) ' = [];']);
	end;
end;


if fmriTEST == 0
  % For each subject, recursively call the program to perform the
  % registration.
  %-----------------------------------------------------------------------
  spm('Pointer','Watch')
  for i=1:nsubjects,
	spm('FigName',['Coregister: working on subj ' num2str(i)],Finter,CmdLine);
	fprintf('\rCoregistering Subject %d: ', i,'\n');

	eval(['PGF    =    PGF' num2str(i) ';']);
	eval(['PFF    =    PFF' num2str(i) ';']);
	eval(['others = others' num2str(i) ';']);

	if p == 1 | p == 3,
		spm_coregister(PGF, PFF, PGG, PFG, others,flags);
	end;
	if p == 2 | p == 3,
		% Write the coregistered images
		%-----------------------------------------------------------------------
		P = str2mat(PGF(1,:),PFF);
		if prod(size(others))>0,
			P = str2mat(P,others);
		end;
		spm_reslice(P,struct('mask',0,'mean',0,'hold',1,'which',1));
	end;
  end;
  fprintf('\r%60s%s', ' ',sprintf('\b')*ones(1,60));
  spm('FigName','Coregister: done',Finter,CmdLine);
  spm('Pointer');

end

statusfile = [fmriDIR filesep 'LOG' filesep 'COREGISTER.OK'];
f99_log(statusfile,'');

tijd=spm('time');
f99_log(logfile,['Coregistration Script ended at ' mat2str(tijd)]);
disp('-----------------------------------------------------------------');
disp('');
return

%_______________________________________________________________________

%_______________________________________________________________________
function edit_defs
global sptl_UsMtlInfrmtn sptl_QckCrg

sptl_UsMtlInfrmtn = spm_input(...
	['Use Mutual Information Registration?'],...
	'+1', 'm',...
	['Use Default Registration|'...
	 'Use Mutual Information Registration'],...
	[0 1], any(sptl_UsMtlInfrmtn)+1,'batch',{},'UsMtlInfrmtn');

if ~sptl_UsMtlInfrmtn,
	sptl_QckCrg = spm_input(...
		['Do first coreg step only?'],...
		 '+1', 'm',...
		['Do all three coreg steps|'...
		 'Do first coreg step only'],...
		[0 1], any(sptl_QckCrg)+1,'batch',{},'QckCrg');
end;
return;
%_______________________________________________________________________



