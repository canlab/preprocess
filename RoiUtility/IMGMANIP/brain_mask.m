function brain_mask(inname,outname)
% brain_mask(inname,outname)
%
% Masks an image with a canonical smoothed gray matter mask.

mask = which('scalped_avg152T1_graymatter_smoothed.img');
if isempty(mask), error('Cannot find canonical mask image');, end

V = spm_vol(str2mat(mask,inname));
chk = check_spm_mat(V(1).mat,V(2).mat)

if ~chk
        [P,mask] = reslice_imgs(inname,mask,1);
end


Q = spm_imcalc_ui(str2mat(mask,inname),outname,'(i1>0) .* i2');



return

