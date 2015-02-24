function scn_setup()
    global SCN_CONSTANTS;
    SCN_CONSTANTS = 1;

    global MKLINKS;
    MKLINKS = '/Volumes/SCNAlpha/matthew_scripts/mklinks2';

    global MEDCON;
    MEDCON = '/usr/local/xmedcon/bin/medcon';
    if(~exist(MEDCON, 'file'))
        MEDCON = '/sw/bin/medcon';
    end

    global FSLDIR;
    FSLDIR = '/usr/local/fsl';
    if(~exist([FSLDIR '/bin/fslchfiletype'], 'file'))
        FSLDIR = '/Applications/FSL/fsl';
    end

    global ANA4DTO3D
    ANA4DTO3D = '/sw/bin/ana4dto3d';
    if(~exist(ANA4DTO3D, 'file'))
        ANA4DTO3D = '/opt/local/bin/ana4dto3d';
    end

    global AFNIDIR
    AFNIDIR = '/usr/local/afni';
    
    global defaults
    spm_defaults();
end