% newSPMs = change_spm_path(SPM_files, old_path, new_path, [save_file])
%
% Goes through the fields containing paths of an SPM.mat and changes them
% Useful when you move an SPM to another location.
% Replace SPM_files with [] to prompt to select SPM.mat files
% NB: Does not overwrite old SPM.mat files by default
%
% old_path = '/Volumes/old_root/';
% new_path = '/Volumes/new_root/';
% SPM_files = '/Volumes/new_root/analysis/subj_whatever/SPM.mat';
% newSPM = change_spm_path(SPM_files, old_path, new_path, 1) % to overwrite existing SPM.mat file
% or
% newSPM = change_spm_path(SPM_files, old_path, new_path) % to just return a structure and not overwrite the old SPM.mat

function newSPMs = change_spm_path(SPM_files, old_path, new_path, save_file)
    if(~exist('save_file', 'var') || isempty(save_file))
        save_file = 0;
    end
    if(isempty(SPM_files))
        SPM_files = spm_select(Inf, 'mat', 'Select SPM.mat files...', [], pwd(), '^SPM\.mat$');
    end
    SPM_files = cellstr(SPM_files);
    
    if(~ischar(old_path) || ~ischar(new_path))
        error('Inputs are not strings.');
    end
    
    for i=1:length(SPM_files)
        SPM_file = SPM_files{i};
        load(SPM_file);
        newSPMs(i) = struct_strrep(SPM, old_path, new_path);

        if(save_file)
            SPM = newSPMs(i);
            save(SPM_file, 'SPM');
        end
    end
end