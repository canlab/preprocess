function changeVox(vox_x,vox_y,vox_z)
% changeVox(vox_x,vox_y,vox_z)
%
% vox_x,y, & z are the new voxel dimensions for a set of img header files.
%
% run in main subject dir:
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

dd = dir('scan*');  % look for scan? directories, where ? is an integer
for i = 1:length(dd)
    if dd(i).isdir  % make sure it a directory
        
        disp(['changing vox size in directory ' dd(i).name])
        [PP,dirs] = spm_list_files(dd(i).name,'ow*img'); % image names are ow?_scan?_vol????.img
        [DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = spm_hread([dd(i).name filesep deblank(PP(1,:))]);
        disp(['   Changing original vox size: ' num2str(VOX(1)) ' ' num2str(VOX(2)) ' '  num2str(VOX(3))]);
        disp(['                           to: ' num2str(vox_x) ' ' num2str(vox_y) ' '  num2str(vox_z)]);
        for j = 1:size(PP,1)
            sts = spm_header_edit('set',[dd(i).name filesep deblank(PP(j,:))],['VOX=[' num2str(vox_x) ' ' num2str(vox_y) ' ' num2str(vox_z) '];']);
        end
        
    end
end

return



