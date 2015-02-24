function use_spm(spm_ver)
    clear global defaults
    clear global SPM*
    clear global st
    clear global transv
    clear global Users
    clear global ChannelFlag
    clear global Verbose
    clear global UFp

    switch(lower(spm_ver))
        case {'5' 'spm5'}
            strip_path_dirs('spm99');
            strip_path_dirs('spm2');
            add_spm_to_path('spm5');
            spm('defaults', 'FMRI');
        case {'2' 'spm2'}
            strip_path_dirs('spm99');
            strip_path_dirs('spm5');
            add_spm_to_path('spm2');
            spm_defaults();
        case {'99' 'spm99'}
            strip_path_dirs('spm2');
            strip_path_dirs('spm5');
            add_spm_to_path('spm99');
    end
    strip_svn_dirs();
end


function add_spm_to_path(spm_ver)
    SPM99_HOME = '/Volumes/SCNAlpha/matlab_code_external/spm99';
    SPM2_HOME = '/Volumes/SCNAlpha/matlab_code_external/spm2';
    SPM5_HOME = '/Volumes/SCNAlpha/matlab_code_external/spm5';

    pathdirs = explode(path(), pathsep);
    wh_spm = find(wh_str(['/' spm_ver], pathdirs), 1);
    if(~isempty(wh_spm))
        addpath(genpath(pathdirs{wh_spm}));
    else
        switch(spm_ver)
            case 'spm5'
                addpath(genpath(SPM5_HOME));
            case 'spm2'
                addpath(genpath(SPM2_HOME));
            case 'spm99'
                addpath(genpath(SPM99_HOME));
        end
    end
end