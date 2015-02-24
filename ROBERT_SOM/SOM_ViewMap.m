
%PAnatomy = spm_get([0 1],'*.img','Pick the anatomy');
PAnatomy = which('scalped_single_subj_T1.img');,
if length(PAnatomy) < 0
    return
end

PSOM     = spm_get([0 1],'*.img','Pick the SOM Image');

if length(PSOM) < 0
    return
end

spm_check_registration(PAnatomy);

mycap = colormap([gray(64);jet(64)]);

spm_orthviews('addtruecolourimage',1,PSOM,mycap,.4,1,0);
