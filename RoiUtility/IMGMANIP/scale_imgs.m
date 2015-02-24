function Pout = scale_imgs(hP,divP)
% Pout = scale_imgs(hP,divP)
%
% Takes a string matrix of image file names hP
% divides each voxel by values in divP
%
% Writes SC* images (SCaled)
%
% tor wager

V = spm_vol(hP);
V2 = spm_vol(divP);

fprintf(1,'\nloading imgs.')
v = spm_read_vols(V);
fprintf(1,'loading divisor.')
v2 = spm_read_vols(V2);

v2 = repmat(v2,[1 1 1 size(v,4)]);

fprintf(1,'scaling\n')
warning off
v2(v2==0) = NaN;
vout = v ./ v2;
warning on

for i = 1:size(hP,1)

	[d,f,e] = fileparts(V(i).fname);
	V(i).fname = [d filesep 'SC' f e];

	vtmp = squeeze(vout(:,:,:,i));

	if i == 1, Pout = V(i).fname;, else, Pout = str2mat(Pout,V(i).fname);, end

	spm_write_vol(V(i),vtmp);
	disp(['Written ' V(i).fname])

end
