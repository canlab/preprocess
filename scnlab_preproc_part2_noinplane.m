function scnlab_preproc_part2_noinplane(varargin)
%scnlab_preproc_part2_noinplane(varargin)
%
% Series of steps to take images from realigned (ravols or MCavols) to
% smoothed, normalized functional images 
%
%
% obj = 'structural/T1.img';
% func = 'r*/ravol*img';
% nfunc = 'r*/wravol*img';
% nobj = 'structural/wT1.img';
% templ = which('avg152T1.img');
% 
% %for multiple subjects, e.g., 
% subjdirs = {'rea28' 'rea29' 'rea30'};
% for i=1:length(subjdirs) % e.g.
%    cd(subjdirs{i});
%    scnlab_preproc_part2_noinplane('obj',obj,'func',func,'nobj',nobj,'nfunc',nfunc,'templ',templ);
%    cd('..');
% end
% 
%

obj = 'structural/T1.img';
func = 'r*/ravol*img';
nfunc = 'r*/wravol*img';
nobj = 'structural/wT1.img';
templ = which('avg152T1.img');

%process var args
if length(varargin) > 0

    for i = 1:length(varargin)
        if strcmp(varargin{i},'targ'), targ = varargin{i+1};,end
        if strcmp(varargin{i},'obj'), obj = varargin{i+1};,end
        if strcmp(varargin{i},'func'), func = varargin{i+1};,end
        if strcmp(varargin{i},'nobj'), nobj = varargin{i+1};,end
        if strcmp(varargin{i},'nfunc'), nfunc = varargin{i+1};,end
        if strcmp(varargin{i},'templ'), templ = varargin{i+1};,end
    end

end


% set the origin of the hi-res SPGR
set_hdr_current_coords(obj);

% get list of functional files
p = get_filename2(func);

% set the origin of the first functional and adjust all functionals
set_hdr_current_coords(p(1,:),p);


spm_check_registration(str2mat(p(1,:),obj))
disp('Check whether images roughly match!!');
s = input('Press a key to go on, or ctrl+c to break');

% coregister in-plane T1 to EPI refVol
% scnlab_coreg_anat2funct(p(1,:),targ);   
scnlab_coreg_anat2funct(p(1,:),obj);   

% old
% % coregister the SPGR to the in-plane T1
% scnlab_coreg_spgr2inplane(targ,obj);


% check the results
spm_check_registration(str2mat(obj,p(1,:)))
spm_print


% normalize to template, apply to functionals, and smooth
pp = scnlab_spm2_norm(1,1,1,1,obj,templ,p,nfunc);


% check results
%template_file = which('avg152t1.img'); % which('T1.mnc'); %%'scalped_avg152T1.img'); %old
spm_check_registration(str2mat(templ,nobj,pp(1,:)))
spm_print

return


