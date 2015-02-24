function adj_mask_coreg(P,p)
% adj_mask_coreg(anat_img_name,funct_img_name)
%
%    1st arg: anat_img_name points to the anatomical image,
%    2nd arg: funct_img_name points to the 1st functional image (or volume)
%
% Return values:
%               None.
%
% This routine will allow you to spatially adjust the position of a structural (anatomical) image
% before co-registering it to it's cooresponding first functional image.
% 
% Some variables:
%
% P = brain extracted hi-res T1
% p = the reference volume, the very first functional of the experiment
%
% the final coregistered anatomical will be m[P], e.g., mesT1_s5.img, if you
% run the full version.  If you enter a 3rd argument, you will run the
% 'expedited' version, and the final output will be P, esT1_s5.img.
%
% This function coregisters without reslicing, to preserve the spatial
% definition of the anatomical.
%
%$ * Project: fMRI Imaging lab: Wager/Ochsner *
%$ * Author: Christopher Currie	(modified from scnlab_coreg_anat2funct.m by Wager)*
%$ * Columbia University	* 
%$ * Department of Psychology FAX: +1 212 854-3609 * 
%$ * 1190 Amsterdam Avenue, Schermerhorn Hall, Room 406 TEL +1 212 854-3608	*
%$ * New York, NY 10003 EMail: ccurrie@paradox.psych.columbia.edu	*
%$ * *
%$ * Version 1.0, March 31 2004.

mytarget = p;
myobject = P;

% s = input('Press return to continue');
%
% simply check the registration
%
spm_check_registration(str2mat(mytarget,myobject))
%
% Display iterative commands to user
%
s = input('Press return to continue');
disp(['YOU ARE IN DIRECTORY: ' pwd])
disp('Use display to shift anatomical esT1_s5.img until the two match')
disp('TYPE: spm_image(''init'',myobject)')
disp('TYPE: spm_check_registration(str2mat(mytarget,myobject))')
disp('TYPE: spm_coregister(mytarget,myobject)')
disp('Type return when they match')
%
% KEYBOARD, when placed in an M-file, stops execution of the file
% and gives control to the user's keyboard. 
% The keyboard mode is terminated by executing the command RETURN
%
keyboard

%
% give 'smooth_and_mask' the hi-res structural image and 
% assign it to Q.
% The routine smooths the anatomical and tries to put a mask
% around what it thinks is the brain in order to remove
% artifacts, such as the skull.  What's left is "just" the anatomy.
%
Q = smooth_and_mask(P);
%
% Two step clean-up: 
%      1) Q (the mask, derived from the anatomical) is applied 
%      to the functional.
%      2) turns zeros into NaN values in the functional image
%      NaN is ignored by the co-registration process, and 
%      function should contain no zero values.
%
out = nanmask(p,Q);
%
% the new functional (Q) is assigned to 'mytarget'
%
P2 = nanmask(P);
mytarget = out;
P2 = nanmask(P);
myobject = P2;

mi_coreg_plugin3



return
