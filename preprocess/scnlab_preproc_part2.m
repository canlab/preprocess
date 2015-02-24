%scnlab_preproc_part2(varargin)
%
% Series of steps to take images from realigned (ravols or MCavols) to
% smoothed, normalized functional images
%
%
% obj = 'structural/T1.img';
% targ = 'structural/T1inplane.img';
% func = 'r*/ravol*img';
% nfunc = 'r*/wravol*img';
% nobj = 'structural/wT1.img';
% templ = which('avg152T1.img');
%
% for i=1:length(subjdirs)
%    cd(subjdirs{i});
%    scnlab_preproc_part2('targ', targ, 'obj', obj, 'func', func, 'nobj', nobj, 'nfunc', nfunc, 'templ', templ);
%    cd('..');
% end

function scnlab_preproc_part2(varargin)
    writefunc = 1;
    dosmooth = 1;
    obj = 'structural/T1.img';
    targ = 'structural/T1inplane.img';
    func = 'r*/ravol*img';
    nfunc = 'r*/wravol*img';
    nobj = 'structural/wT1.img';
    templ = which('avg152T1.img');

    %process var args
    if length(varargin) > 0
        for i = 1:length(varargin)
            if ischar(varargin{i})
                if strcmp(varargin{i}, 'targ'), targ = varargin{i+1}; varargin{i+1} = []; end
                if strcmp(varargin{i}, 'obj'), obj = varargin{i+1}; varargin{i+1} = []; end
                if strcmp(varargin{i}, 'func'), func = varargin{i+1}; varargin{i+1} = []; end
                if strcmp(varargin{i}, 'nobj'), nobj = varargin{i+1}; varargin{i+1} = []; end
                if strcmp(varargin{i}, 'nfunc'), nfunc = varargin{i+1}; varargin{i+1} = []; end
                if strcmp(varargin{i}, 'templ'), templ = varargin{i+1}; varargin{i+1} = []; end

                if strcmp(varargin{i}, 'nosmooth'), dosmooth = 0; end
                if strcmp(varargin{i}, 'nofunc'), writefunc = 0; end
                if strcmp(varargin{i}, 'funclist'), p = varargin{i+1}; varargin{i+1} = []; end
            end
        end
    end


    % set the origin of the hi-res SPGR
    set_hdr_current_coords(obj);

    % get list of functional files if we don't have it already
    if exist('p', 'var')
        % skip getting filenames; we have them
    else
        p = get_filename2(func);
    end

    % set the origin of the in-plane T1 and adjust all functionals
    set_hdr_current_coords(targ, p);

    spm_check_registration(str2mat(targ, obj, p(1,:)))
    disp('Check whether images roughly match!!');
    s = input('Press a key to go on, or ctrl+c to break');

    % coregister in-plane T1 to EPI refVol
    scnlab_coreg_anat2funct(p(1,:), targ);

    % coregister the SPGR to the in-plane T1
    scnlab_coreg_spgr2inplane(targ, obj);


    % check the results
    spm_check_registration(str2mat(obj, targ, p(1,:)))
    spm_print


    % normalize to template, apply to functionals, and smooth
    pp = scnlab_spm2_norm(1, 1, writefunc, dosmooth, obj, templ, p, nfunc);


    % check results
    if(isempty(pp))
        spm_check_registration(str2mat(templ, nobj))
    else
        spm_check_registration(str2mat(templ, nobj, pp(1,:)))
    end
    spm_print
end