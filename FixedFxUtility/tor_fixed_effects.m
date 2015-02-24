function [vol,maskingImage,V] = tor_fixed_effects(file_loc,crit_p,cl_size)
% [vol,maskingImage,V] = tor_fixed_effects(file_loc,crit_p,cl_size)
% shell for running fixedfx2, fixed effects analysis
%
% file_loc is the location of files, including image name
% as you'd enter if you were listing them with ls.
%
% example:
% '/data/biman4/RESULTS/model2/sub*/spmT_0002.img'

% ------------------------------------------------------
% * get files
% ------------------------------------------------------

% Windows
fnames = getfiles2(file_loc);
nimgs = size(fnames,1);

% Unix / Linux
%[fnames,nimgs] = getfiles(file_loc);
%fnames = str2mat(fnames);

if isempty(fnames), error('Couldn''t find any files.'),end

disp(['Found ' num2str(nimgs) ' images for fixed effects analysis.'])

fnames

% ------------------------------------------------------
% * do fixed effects analysis
% ------------------------------------------------------

[vol,maskingImage,V] = fixedfx2(fnames,cl_size,crit_p);


% ------------------------------------------------------
% * display results
% ------------------------------------------------------

str = input('Display results? (y\\n) ','s');

if strcmp(str,'y')

   if ~isempty(vol)
	try
		disp('Choose canonical image.')
		spm_image
		spm_orthviews('AddBlobs',1,V.XYZ,V.Z,V.mat)

		%disp('Click on white space for text')
		%ginput(1)

%pos = [2.5556    6.5314
%    2.4815    5.5745
%    2.5556    4.7279
%    2.4815    3.7341
%    2.4815    3.0348
%    2.5556    2.2987
%    2.3333    1.5626
%   2.4815    0.6792];
		disp('Click 8 white spaces for text lines')
		pos = ginput(8);

		text(pos(1,1),pos(1,2),V.descrip(end-40:end))
		text(pos(2,1),pos(2,2),'fixed effects analysis')
		text(pos(3,1),pos(3,2),['p < ' num2str(V.crit_p)])
		text(pos(4,1),pos(4,2),['t > ' num2str(V.crit_t)])
		text(pos(5,1),pos(5,2),['cluster size >= ' num2str(V.cl_size)])
		text(pos(6,1),pos(6,2),['Clusters: ' num2str(V.numClusters)])
		text(pos(7,1),pos(7,2),V.fname(1:end-20))
		text(pos(8,1),pos(8,2),V.fname(end-20:end))
		
	catch
		disp('display problem')
	end

   end
end

	% ------------------------------------------------------
	% * get clusters for each image
	% ------------------------------------------------------
	
str = input('Extract clusters? (y\\n) ','s');

if strcmp(str,'y')
	
		fn2{1} = fnames;
		%fn2
		VOL.VOX = V.voxSize;
		[clusters] = tor_extract_rois(fn2,V,VOL);
		% save V in mat file
		% ------------------------------------------------------
		[p,f,e] = fileparts(V.fname);
		f = ['V_' f];
		fname = [p filesep f];
		disp(['Saving V information and clusters in: ' fname '.mat'])
		eval(['save ' fname ' V clusters'])
	
		%disp('Problem getting clusters: skipping.')

else
	a = V.fname;
	[b,a,ext] = fileparts(a);
	a = [a ext];
	clusters(1).title = a;
end

	% ------------------------------------------------------
	% * write results mask for fixed effects analysis
	% ------------------------------------------------------

	V2 = V;
	V2.fname = ['res_mask_' clusters(1).title];
	disp(['Writing binary mask to ' V2.fname])
	V = spm_write_vol(V2,V.maskingImage);

end

return