function scnlab_coreg_anat2funct(P,p,varargin)
% scnlab_coreg_anat2funct(anat_img_name,funct_img_name,[suppress iterative])
%
% P = brain extracted hi-res T1
% P ='/Users/scnlab/Kosslyn/Data_and_Tools/IMAGING_DATA/Amygdala_Face_Class_Data/ow3/anatomy/esT1_s5.img'
% p = the reference volume, the very first functional of the experiment
% p ='/Users/scnlab/Kosslyn/Data_and_Tools/IMAGING_DATA/Amygdala_Face_Class_Data/ow3/scan1/firstTR/refVol.img'
%
% the final coregistered anatomical will be m[P], e.g., mesT1_s5.img, if you
% run the full version.  If you enter a 3rd argument, you will run the
% 'expedited' version, and the final output will be P, esT1_s5.img.
%
% This function coregisters without reslicing, to preserve the spatial
% definition of the anatomical.

mytarget = p;
myobject = P;
%mi_coreg_plugin3

if length(varargin) > 0
    
else
    s = input('Press return to continue');
    spm_check_registration(str2mat(mytarget,myobject))
    s = input('Press return to continue');
    pwd
    disp('Use display to shift anatomical inplane_T1.img until the two match')
    disp('TYPE: spm_image(''init'',myobject)')
    disp('TYPE: spm_check_registration(str2mat(mytarget,myobject))')
    disp('TYPE: spm_coregister(mytarget,myobject)')
    disp('Type return when they match')

    keyboard

    Q = smooth_and_mask(P);
    out = nanmask(p,Q);

    P2 = nanmask(P);
    mytarget = out;
    P2 = nanmask(P);
    myobject = P2;

    mi_coreg_plugin3
end


return
