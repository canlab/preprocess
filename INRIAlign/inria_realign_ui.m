function inria_realign_ui(arg1)
% User Interface for inria_realign.
%___________________________________________________________________________
%
% The INRIAlign toolbox enhances the standard SPM realignment routine
% (see topic: spm_realign_ui). In the latter, rigid registration is
% achieved by minimization of the sum of squared intensity differences
% (SSD) between two images. As noted by several SPM users, SSD based
% registration may be biased by a variety of image artifacts and also
% by activated areas. To get around this problem, INRIAlign reduces
% the influence of large intensity differences by weighting errors
% using a non-quadratic, slowly-increasing function (rho
% function). This is basically the principle of an M-estimator.
%
% When launching INRIAlign, the user may select a specific rho
% function as well as an associated relative cut-off distance (which
% is needed by most of the rho functions). By default, the rho
% function is that of Geman-McClure while the relative cut-off
% distance is set to 2.5.
%
% Apart from this distinction, the method is very similar to
% spm_realign and uses the same editable default parameters. Most of
% the implementation has been directly adapted from the code written
% by J. Ashburner (See spm_realign).
%
% This version is updated for SPM2 and is no longer SPM99 compatible. 
% Note: it does not work in Batch mode. 
%
%__________________________________________________________________________
% Refs:
%
% Freire L, Roche A & Mangin JF (2002). What is the best similarity 
% measure for motion correction in fMRI time series? 
% IEEE Trans. Med. Imag., 21(5):470--484.
%
% Rousseeuw PJ and Leroy AM (1987). Robust Regression and Outlier
% Detection. Wiley Series in Probability and Mathematical
% Statistics. 
%
%__________________________________________________________________________
% @(#)inria_realign_ui.m  2.0  Alexis Roche (INRIA/CEA, France) 05/11/23

global defaults

inria_defaults = edit_local_defaults( 'Defaults' );

% User interface.
%_______________________________________________________________________
SPMid = spm('FnBanner',mfilename,'1.1');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','INRIAlign');
spm_help('!ContextHelp','inria_realign_ui.m');

pos = 1;

% Enables editing default parameters
%-----------------------------------------------------------------------
if spm_input(['Check local parameters?'], pos, 'm',['Use defaults|' ...
	      'Edit parameters'], [0 1], 1);
  pos = pos + 4;
  inria_defaults = edit_local_defaults( 'User' );

end,
	
n = spm_input('number of subjects', pos, 'e', 1 ); 
if (n < 1)
	spm_figure('Clear','Interactive');
	return;
end

P = cell(n,1);
pos = pos + 1;
for i = 1:n,
	if strcmp(lower(defaults.modality), 'fmri'),
		ns = spm_input(['num sessions for subject ',num2str(i)], pos, 'e', 1); 
		pp = cell(1,ns);
		for s=1:ns,
			p = '';
			while size(p,1)<1,
				p = spm_get(Inf,'.img',...
				['scans for subj ' num2str(i) ', sess' num2str(s)]);
			end;
			pp{s} = p;
		end;
		P{i} = pp;
	else, %- no batch mode for 'PET'
		p  = cell(1,1);
		p{1} = '';
		while size(p{1},1)<1,
		      p{1} = spm_get(Inf,'.img',...
			  ['select scans for subject ' num2str(i)]);
		end;
		P{i} = p;
	end;
end;

% Instantiate parameter struct for inria_realign.m
FlagsC = struct( 'quality', defaults.realign.estimate.quality, 'fwhm', inria_defaults.Smooth, ...
                 'rho_func', inria_defaults.CostFun, 'cutoff', inria_defaults.CutOff );
%% Registration to mean allowed for PET
if strcmp(lower(defaults.modality),'pet'),
  FlagsC.rtm = []; 
end;

FlagsC,

WhchPtn = spm_input('Which option?', pos, 'm',...
	            'Coregister only|Reslice Only|Coregister & Reslice',...
	            [1 2 3], 3); 
pos = pos + 1;

PW = '';
if (WhchPtn == 1 | WhchPtn == 3) & defaults.realign.estimate.weight,
	if spm_input(...
		['Weight the reference image(s)?'],...
		2, 'm',...
		['Dont weight registration|'...
		 'Weight registration'], [0 1], 1);
		PW = spm_get(n,'.img',...
			'Weight images for each subj');
	end;
end;

% Reslicing options
%-----------------------------------------------------------------------
if WhchPtn == 2 | WhchPtn == 3,
	
	FlagsR = struct('interp',defaults.realign.write.interp,...
		        'wrap',defaults.realign.write.wrap,...
		        'mask',defaults.realign.write.mask,...
		        'which',2,'mean',1);
	if strcmp(lower(defaults.modality),'pet'), FlagsR.wrap = [0 0 0]; end;

FlagsR, 

	p = spm_input('Create what?','+1','m',...
		[' All Images (1..n)| Images 2..n|'...
		 ' All Images + Mean Image| Mean Image Only'],...
		[1 2 3 4], 3);
	if p==1, FlagsR.which = 2; FlagsR.mean = 0; end
	if p==2, FlagsR.which = 1; FlagsR.mean = 0; end
	if p==3, FlagsR.which = 2; FlagsR.mean = 1; end
	if p==4, FlagsR.which = 0; FlagsR.mean = 1; end
	
end;	
	
% Do the job
spm('Pointer','Watch');
for i = 1:n
	spm('FigName',['INRIAlign: working on subject ' num2str(i)],Finter,CmdLine);
	fprintf('\rRealigning Subject %d: ', i);
	if WhchPtn==1 | WhchPtn==3,
		if ~isempty(PW), FlagsC.PW = deblank(PW(i,:)); end;
		inria_realign(P{i}, FlagsC);
	end
	if WhchPtn==2 | WhchPtn==3,
		spm_reslice(P{i},FlagsR)
	end;
end
fprintf('\r%60s%s', ' ',sprintf('\b')*ones(1,60));
spm('FigName','INRIAlign: done',Finter,CmdLine);
spm('Pointer');

return;


%_______________________________________________________________________
%_______________________________________________________________________
function LDefs = edit_local_defaults( arg ) 

defSmooth = 0; 
defCostFun = 'geman'; 
defCutOff = 2.5; 

LDefs = struct( 'Smooth', defSmooth, 'CostFun', defCostFun, 'CutOff', defCutOff ); 

if strcmp( arg, 'User' ), 

  LDefs.Smooth = spm_input('Spatial smoothing (mm)', '+1', 'e', defSmooth); 

  tmp2 = str2mat('quadratic','absolute','huber','cauchy','geman','leclerc','tukey');
  tmp = strmatch(defCostFun,tmp2);
  aux = spm_input('Cost function?','+1','m',...
	 	  ['Quadratic  (SPM standard)|Absolute value|Huber|'...
		   'Cauchy|Geman-McClure|Leclerc-Welsch|Tukey|'], (1:size(tmp2,1)), tmp ); 

  LDefs.CostFun = deblank(tmp2(aux,:));

  if strcmp(LDefs.CostFun,'quadratic') | strcmp(LDefs.CostFun,'absolute'), 
    LDefs.CutOff = Inf; 
    return; 
  end,

  LDefs.CutOff = spm_input('Relative cut-off distance', '+1', 'e', defCutOff); 

end, 

return;

