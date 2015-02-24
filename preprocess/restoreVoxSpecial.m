function restoreVoxSpecial
% 
% Use only with ow0 through ow6 to restore original vox size in y-dimension
% 
% To insure that the anatomical and functionals are NOT flipped along the
% x-dim. relative to each other, 
% 1st) run spm_defaults.m after starting matlab, i.e., type "spm_defaults"
% in command window.
% 2nd) start spm ("spm fmri")
% 3rd) type "defaults.analyze.flip" -- SHOULD BE SET TO 0!
% 4th) Both the refVol (and all functionals) and the anatomical should have
%      the same size for the x voxel dimension, i.e.,
% [DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = spm_hread(functional) &
% [DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = spm_hread(anatomy);
%      should return the same SIGN (i.e., positive) for each X, Y, Z
%      components of VOX.
%
%$ * Project: fMRI Imaging lab: Wager/Ochsner *
%$ * Author: Christopher Currie	*
%$ * Columbia University	* 
%$ * Department of Psychology FAX: +1 212 854-3609 * 
%$ * 1190 Amsterdam Avenue, Schermerhorn Hall, Room 406 TEL +1 212 854-3608	*
%$ * New York, NY 10003 EMail: ccurrie@paradox.psych.columbia.edu	*
%$ * *
%$ * Version 1.0, April 1 2004.

firstTR = 'scan1/firstTR/refVol.img';
%
sts = spm_header_edit('set',firstTR,['VOX=[2.9688 2.9688 4]']);
if sts == 0
   disp(['sts = ' num2str(sts) ', spm_header_edit was successful for ' firstTR]);
else
   disp(['sts = ' num2str(sts) ', spm_header_edit was NOT successful for ' firstTR ]);
end

[DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = spm_hread('scan1/firstTR/refVol.img');
disp(['VOX = ' num2str(VOX)]);

dd = dir('scan*');
for i = 1:length(dd)
    if dd(i).isdir
        
        disp(['Un-shrinking ' dd(i).name])
        [PP,dirs] = spm_list_files(dd(i).name,'ow*img');
        
        for j = 1:size(PP,1)
            sts = spm_header_edit('set',[dd(i).name filesep deblank(PP(j,:))],['VOX=[2.9688 2.9688 4];']);
        end
        
    end
end

return



