%preproc_part2_spm2(varargin)
%
% Series of steps to take images from realigned (ravols) to
% smoothed, normalized functional images (swravols)
%
% Keyword Parameters:
%    anat - high-resolution, T1 anatomical image
%    func - path string for filenames() for the functional images
%
% Optional Keyword Parameters:
%    inplane - T2* image with the same number of slices as a functional image, 
%        but much higher resolution within the slice - useful as an intermediate image
%    templ - caconical image to warp to - (default: avg152T1.img);
% 
% Flag Parameters (not followed by any values)
%    nofunc - turns off normalization of functional images
%    nosmooth - turns off smoothing of functional images
%
% E.g.:
%
% anat = 'structural/T1.img';
% func = 'run*/ravol*img';
%
% %for multiple subjects, e.g.,
% subjdirs = {'rea28' 'rea29' 'rea30'};
% for i=1:length(subjdirs) 
%    cd(subjdirs{i});
%    preproc_part2_spm2('anat', anat, 'func', func);
%    cd('..');
% end
%
% % with an inplane .img:
% inplane = 'structural/T1inplane.img';
% preproc_part2_spm2('anat', anat, 'inplane', inplane, 'func', func);
%
% % without smoothing:
% preproc_part2_spm2('anat', anat, 'func', func, 'nosmooth');

function preproc_part2_spm2(varargin)
    use_spm('SPM2');
    writefunc = 1;
    dosmooth = 1;
    anat = 'structural/T1.img';
    inplane = [];
    func = 'r*/ravol*img';
    templ = which('avg152T1.img');

    %process var args
    if ~isempty(varargin)
        for i = 1:length(varargin)
            if(ischar(varargin{i}))
                switch(varargin{i})
                    case {'inplane' 'targ'}
                        inplane = varargin{i+1};
                    case {'anat', 'anat_img', 'anat'}
                        anat = varargin{i+1};
                    case 'func'
                        func = varargin{i+1};
                    case 'templ'
                        templ = varargin{i+1};
                    case 'nosmooth'
                        dosmooth = 0;
                    case 'nofunc'
                        writefunc = 0;
                end
            end
        end
    end

    %if(iscell(anat) && iscell(func) && (iscell(func{1}) || ischar(func{1}))

    % set the origin of the hi-res SPGR
    set_hdr_current_coords(anat);

    % get list of functional files
    funcs = filenames(func);
    [pathstr, name, ext] = fileparts(func);
    wfunc = fullfile(pathstr, ['w' name ext]);

    % set the origin of the inplane if it exists, otherwise the first functional, and adjust all functionals
    if(~isempty(inplane))
        set_hdr_current_coords(inplane, char(funcs));
    else
        set_hdr_current_coords(funcs{1}, char(funcs));
    end

    spm_check_registration(strvcat(anat, inplane, funcs{1}))  %#ok
    disp('Check whether images roughly match!!');
    s = input('Press a key to go on, or ctrl+c to break');

    if(isempty(inplane))
        % coregister hires T1 to the reference func volume
        scnlab_coreg_anat2funct(funcs{1}, anat);
    else
        % coregister in-plane T1 to the reference func volume
        % and then coregister the hires T1 to the in-plane T1
        scnlab_coreg_anat2funct(funcs{1}, inplane);
        scnlab_coreg_spgr2inplane(inplane, anat);
    end

    % check the results
    spm_check_registration(strvcat(anat, inplane, funcs{1}));
    spm_print();

    % normalize to template, apply to functionals, and smooth
    normalized_files = scnlab_spm2_norm(1, 1, writefunc, dosmooth, anat, templ, char(funcs), wfunc);

    % check results
    [pathstr, name, ext] = fileparts(anat);
    wanat = fullfile(pathstr, ['w' name ext]);
    spm_check_registration(str2mat(templ, wanat, normalized_files(1,:)));
    spm_print();
end