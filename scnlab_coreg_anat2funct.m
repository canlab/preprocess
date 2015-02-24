function scnlab_coreg_anat2funct(targ,obj,varargin)
% scnlab_coreg_anat2funct(funct_img_name(target),anat_img_name(obj),[suppress iterative])
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

global defaults
spm_defaults

coreg(targ,obj)

if length(varargin) > 0
    
else
    s = 0;
    while ~(strcmp(s,'done'))
        s = input('Type done if finished, or press return to adjust anatomical: ','s');
        
        if ~strcmp(s,'done')
            pwd
            disp('Use display to shift anatomical inplane_T1.img until the two match')
            disp(['TYPE: spm_image(''init'', obj)'])
            disp('TYPE: spm_check_registration(str2mat(obj, targ))')
            %disp('spm99 -- TYPE: spm_coregister(mytarget,myobject)')
            disp('Type return when they match')

            keyboard

            coreg(targ,obj)
            
        end
    end
end

return






function coreg(targ,obj)

global defaults

x  = spm_coreg(spm_vol(targ), spm_vol(obj));

PO = [];    % other images

		if isempty(PO),
			PO = obj;
		else,
			PO = str2mat(obj,PO);
		end;

%write out mat files; also, other images
M  = inv(spm_matrix(x));
		MM = zeros(4,4,size(PO,1));
		for j=1:size(PO,1),
			MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
		end;
		for j=1:size(PO,1),
			spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
		end;


spm_check_registration(str2mat(obj,targ))


return
