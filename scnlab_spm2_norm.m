% pp = scnlab_spm2_norm(donorm, writeanat, writefunc, dosmooth, obj, templ, applytofunclist, [norm func wildcard], [sample at functional resolution?: default 1])
%
% donorm:   1 or 0, do normalization of anatomical and write
% writeanat: write normalized anatomical/structural image
% writefunc: write normalized functional images
% dosmooth: smooth all functional images
%
% templ is the template to normalize to
% obj is name of hi-res SPGR to apply to
% applyto is names of all functionals to apply norm to
%
% Example:
% x = scnlab_spm2_norm(1, 1, 1, 1, obj, p);

function pp = scnlab_spm2_norm(donorm, writeanat, writefunc, dosmooth, obj, templ, applyto, varargin)

    global defaults
    spm_defaults();
    
    nfunc = 'scan*/w*vol*img';
    resampling_at_func_res = 1;

    if length(varargin) > 0
        nfunc = varargin{1};
    end
    if length(varargin) > 1
        resampling_at_func_res = varargin{2};
    end


    Vtempl = spm_vol(templ);

    Vobj = spm_vol(obj);
    objparams = spm_imatrix(Vobj.mat);
    objvoxsize = objparams(7:9);

    fwhm = [8 8 8];                 % smoothing kernel
    defaultnormvoxsize = defaults.normalise.write.vox;
    matname = [spm_str_manip(obj, 'sd') '_sn.mat'];


    % anatomical normalization
    % ---------------------------------------------
    if donorm
        % Determine parameters
        spm_normalise(Vtempl, Vobj, matname, [], [], defaults.normalise.estimate);
    end

    % write anatomical
    % ---------------------------------------------
    if writeanat
        defaults.normalise.write.vox = objvoxsize;
        spm_write_sn(Vobj, matname, defaults.normalise.write);
        defaults.normalise.write.vox = defaultnormvoxsize;
    end


    % write functionals
    % ---------------------------------------------
    if writefunc
        if(resampling_at_func_res)
            Vfunc = spm_vol(applyto(1, :));
            funcparams = spm_imatrix(Vfunc.mat);
            funcvoxsize = funcparams(7:9);
            defaults.normalise.write.vox = funcvoxsize;
        end
        
        spm_write_sn(applyto, matname, defaults.normalise.write);
        defaults.normalise.write.vox = defaultnormvoxsize;
    end

    % get normalized file names (w* in spm2, n* in spm 99)
    pp = filenames(nfunc, 'char');
    
    if dosmooth
        % get handles for interactive window
        [Finter, Fgraph, CmdLine] = spm('FnUIsetup', 'Smooth');

        disp('Smoothing normalized functional files')

        % do the smoothing
        % ---------------------------------------------
        spm('Pointer', 'Watch');
        spm('FigName', 'Smooth: working', Finter, CmdLine);
        spm_progress_bar('Init', size(pp, 1), 'Smoothing', 'Volumes Complete');
        for j = 1:size(pp, 1)
            current_file = deblank(pp(j, :));                                % input image
            [path, name, ext, ver] = fileparts(current_file);                  % old path - same as input

            smoothed_name = fullfile(path, ['s' name ext ver]);      % output image
            spm_smooth(current_file, smoothed_name, fwhm);
            spm_progress_bar('Set', j);
        end
        spm_progress_bar('Clear', j);
        spm('FigName', 'Smooth: done', Finter, CmdLine);
        spm('Pointer');
    end
end
