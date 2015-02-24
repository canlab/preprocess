% ------------------------------------------------
% set up window stuff
% ------------------------------------------------

% SPMid = spm('FnBanner',mfilename,'2.8');
% [Finter,Fgraph,CmdLine] = spm('FnUIsetup','MI Coregister');
% spm_help('!ContextHelp','spm_mireg_ui.m');

% spm('FigName',['MI Coreg: working on ' OPT.subjcode],Finter,CmdLine);
% fprintf('\nCoregistering Subject %s\n', OPT.subjcode);

% ------------------------------------------------
% get inputs to spm_mireg_ui
% ------------------------------------------------
mireg = struct('VG',[],'VF',[],'PO','');
		
% select target(s)
%PG = spm_get(1,'.img', ['select target image for subject ' num2str(i)]);
PG = deblank(mytarget);
mireg.VG = spm_vol(PG);
		
% select object(s)
% PF = spm_get(1,'.img', ['select object image for subject ' num2str(i)]);
VF = deblank(myobject);
mireg.VF = spm_vol(VF);

% ------------------------------------------------
% mutual information coregistration
% ------------------------------------------------

x = spm_mireg(mireg.VG, mireg.VF);
M = inv(spm_matrix(x));


% ------------------------------------------------
% print output - 2 pages
% ------------------------------------------------
% spm_print;
% spm_coregister(mytarget, myobject)
% spm_print