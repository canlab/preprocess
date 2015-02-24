function scnlab_coreg_spgr2inplane(targ,obj)

global defaults

spm_defaults

flags.cost_fun = 'ncc'; % w/i modality best

%temp = which('T1.mnc');
%targ = 'anatomy/inplane_T1.img';
%obj = 'anatomy/t1spgr.img';
%spm_coregister(targ,obj,temp,temp);%%%OLD spm99

x  = spm_coreg(spm_vol(targ), spm_vol(obj),flags);

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