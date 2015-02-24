function tor_image_compare_clusters(varargin)
% tor_image_compare_clusters(varargin)
% inputs are cluster sets (vector of cluster structures)

disp('Choose canonical image.')
spm_image

for i = 1:nargin
	color{i} = input(['Enter color vector for contrast ' num2str(i) ': ']);
end

for i = 1:nargin

	clusters = varargin{i};
	XYZ= []; Z = []; XYZsig{1} = []; Zsig{1} = []; XYZsig{2} = []; Zsig{2} = []; 

	fnames = clusters(1).imnames;
	for j = 1:length(clusters)
		XYZ = [XYZ clusters(j).XYZ];
		Z = [Z clusters(j).Z];

		if isfield(clusters(j),'con1_XYZ') & isfield(clusters(j),'con2_XYZ')
			XYZsig{1} = [XYZsig{1} clusters(j).con1_XYZ];
			XYZsig{2} = [XYZsig{2} clusters(j).con2_XYZ];
		end
	end

	Zsig{1} = ones(1,size(XYZsig{1},2));
	Zsig{2} = ones(1,size(XYZsig{2},2));

	disp([clusters(1).title ' :  ' num2str(size(XYZ,2)) ' suprathreshold voxels.'])

	% use the V structure from the 1st input image
	% ------------------------------------------------------
	V = spm_vol(fnames(1,:));

	% ------------------------------------------------------
	% * display results
	% ------------------------------------------------------

	spm_orthviews('AddColouredBlobs',1,XYZ,Z,V.mat,color{i} ./ 3)

	if isfield(clusters(1),'con1_XYZ') & isfield(clusters(1),'con2_XYZ')
		spm_orthviews('AddColouredBlobs',1,XYZsig{i},Zsig{i},V.mat,color{i})
	end

end




return