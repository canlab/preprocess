function Q = just_mask(P)
%Q = just_mask(P)
% masks at 0
%
% used in our proc stream::
% P = image to mask
% Q = the masked image

[path,name,ext]=fileparts(P);
if isempty(path), path = '.'; end
Q = [path filesep name '_mask' ext];
% spm_smooth(P,Q,8);
Q = spm_imcalc_ui(Q,Q,'i1 > .0',{0 0 spm_type('float') 0});

return

