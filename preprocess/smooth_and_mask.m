function Q = smooth_and_mask(P,varargin)
%Q = smooth_and_mask(P,[OPTION: mask and reslice with 2 x 2 x 2 standard brain])
% smooths with 8 mm filter and masks at .5
%
% used in our proc stream::
% P = brain extracted hi-res T1
% P =
% '/Users/scnlab/Kosslyn/Data_and_Tools/IMAGING_DATA/Amygdala_Face_Class_Data/ow3/anatomy/esT1_s5.img'
% p = the reference volume, the very first functional of the experiment
% p =
% '/Users/scnlab/Kosslyn/Data_and_Tools/IMAGING_DATA/Amygdala_Face_Class_Data/ow3/scan1/firstTR/refVol.img'
% Q = smooth_and_mask(P);
% out = nanmask(p,Q);
%
% P2 = nanmask(P)
%mytarget = out
%P2 = nanmask(P)
%myobject = P2
%mi_coreg_plugin3

[d,f,e]=fileparts(P);
if isempty(d), d = '.'; end
Q = [d filesep 'sm_' f e];


if length(varargin) > 0, 
    dostandard = varargin{1};, 
    spm_smooth(P,Q,dostandard);     % use input smoothing kernel
else, dostandard = 0;,
    spm_smooth(P,Q,8);
end

if dostandard
    pmask = which('scalped_avg152T1_graymatter_smoothed.img');
    [oldQ,Q] = reslice_imgs(pmask,Q,0);
    inp = str2mat(Q,pmask);
    Q = spm_imcalc_ui(inp,Q,'i1 .* (i2 > 0)',{0 0 spm_type('int16') 0});

else
    Q = spm_imcalc_ui(Q,Q,'i1 > .5',{0 0 spm_type('float') 0});
end


return

