function sci = sub_cluster_index(cl,sphere)
% sci = sub_cluster_index(cl,sphere)
%
% function: to find an index of which sub-cluster, defined by a local max,
% each voxel in the cluster belongs to (if any)
%
% input: a cluster structure, defined with tor_extract_clusters (e.g.)
%	and a vector of spheres with coordinates, 
%  	from find_sphere by Russ Poldrack
%
% output: an index, which may be appended to the original cluster vector,
%	of the sub-cluster sphere to which each voxel belongs.
%
% this function is used in tor_get_spheres, which breaks a set of 
%	clusters into sub-clusters as defined by spm_max.
%
% tor wager, 12/12/01

sci = zeros(1,size(cl.XYZ,2));
ov_index = 0;

% -----------------------------------------------
% for each voxel, find which sphere it's in
% -----------------------------------------------

for i = 1:size(cl.XYZ,2)

	for j = 1:length(sphere)

		spXYZ = sphere(j).XYZ;

		clXYZ = repmat(cl.XYZ(:,i),1,size(spXYZ,2));

		mytest = sum(spXYZ == clXYZ) == 3;


		if sum(mytest) > 0

			if sci(i) > 0
				ov_index = ov_index + 1;
				%warning(['Voxel ' num2str(cl.XYZ(:,i)') ' index ' num2str(i) ' is in more than one sphere!'])
	 		end

			% -----------------------------------------------
			% record the sphere index in the sci vector
			% -----------------------------------------------
		   	sci(i) = j;

			% -----------------------------------------------
			% warning if something weird is happening
			% -----------------------------------------------

			if sum(mytest) > 1
				warning(['Voxel ' num2str(cl.XYZ(:,i)') ' index ' num2str(i) ' has more than one match in sphere!'])
			end

		end

	end

end

if ov_index > 0
	warning([num2str(ov_index) ' voxels found in more than one sphere! (counting repeats)'])
end

return
