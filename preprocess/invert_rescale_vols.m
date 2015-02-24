function invert_rescale_vols(sf)
% invert_rescale_vols(sf)
%
% This function will restore images (functional & structural) that were 
% reduced along the y-dimension using scnlab_rescale_vols, i.e., just input
% the same scale factor you used with scnlab_rescale_vols.m
%
% start in main subject dir:
% /Users/scnlab/Kosslyn/Data_and_Tools/IMAGING_DATA/Amygdala_Face_Class_Data/ow3
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
% get current voxel dimensions
%
[DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = spm_hread(firstTR);
%
% invert the scaling
%
sts = spm_header_edit('set',firstTR,['VOX=[' num2str(VOX(1)) ' ' num2str(VOX(2)/(1-sf)) ' ' num2str(VOX(3)) ']']);

dd = dir('scan*');
for i = 1:length(dd)
    if dd(i).isdir
        
        disp(['Shrinking ' dd(i).name])
        [PP,dirs] = spm_list_files(dd(i).name,'ow*img');
        
        for j = 1:size(PP,1)
            sts = spm_header_edit('set',[dd(i).name filesep deblank(PP(j,:))],['VOX=[-2.97 ' num2str(2.97*(1-sf)) ' 4];']);
        end
        
    end
end

scnlab_coreg_anat2funct(anatname,P1,1)
return



