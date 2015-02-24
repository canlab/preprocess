function scnlab_rescale_vols(sf,anatname)
% scnlab_rescale_vols(sf,anatname)
%
% scale factor (sf) is (funct-anat) / func, in mm
%
% start in main subject dir:
% /Users/scnlab/Kosslyn/Data_and_Tools/IMAGING_DATA/Amygdala_Face_Class_Data/ow3
%
% sf =~ .96, different for each subject 
%
% ONLY RUN to shrink y-axis if voxels are 2.97 x 2.97 x 4 (ow0-6 ONLY)

if isempty(anatname)
    anatname = 'anatomy/esT1_s5.img';
else
    anatname = ['anatomy/' anatname];
end

P1 = 'scan1/firstTR/refVol.img';
sts = spm_header_edit('set',P1,['VOX=[-2.97 ' num2str(2.97*(1-sf)) ' 4]']); spm_image('init',P1)

spm_check_registration(str2mat(P1,anatname))

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



