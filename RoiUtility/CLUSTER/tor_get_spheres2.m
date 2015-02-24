function [cl_out,clusters,sphere,SPM,VOL] = tor_get_spheres(clusters,radius,SPM,VOL)
%  [cl_out,clusters,sphere,SPM,VOL]] = tor_get_spheres(clusters,radius,SPM,VOL)
%
% function: to create a vector of clusters (cl_out) consisting of thresholded
%	activation clusters broken into sub-clusters based on local maxima.
%	Sub-cluster centers are defined with spm_max, and
%	cluster boundaries are defined by a sphere of a chosen radius, in mm.
%	
%	To be included in a sub-cluster, a voxel must be within the sphere
%	and must be significant at the chosen threshold.
%
% inputs:
%	radius					in mm; if empty, user is asked for it.		
%	SPM, VOL		 from spm_getSPM, stored in SPM.mat file
%							if empty, prompts for them, and thresholds at
%							height read from cluster structure and extent
%							threshold you enter.
%							You can also enter other structures with the
%							required information 
%                           (e.g., V from tor_fixed_effects, CLU structure)
%
%                       required fields:
%                       SPM.XYZ, XYZmm, Z
%                       VOL.VOX, M
%
%	clusters			a vector of cluster structures, created in (e.g.)
%							tor_extract_clusters.m
%
% outputs:
%	cl_out				the new cluster structure with sub-clusters.
%	clusters			the original clusters input, with sub_cluster_index added.
%
% tor_get_spheres2
% 	difference from tor_get_spheres:
%	1st version assumes that clusters were acquired in the same order
%	as spm_clusters will produce on the input SPM structure
%
%	2nd version (2) fixes this limitation.
%	both versions assume that a sphere cannot contain voxels from
%	more than 1 cluster.
%
% by tor wager, 12/13/01, fixed bug 1/30/02 in timeseries reporting
% doesn't work quite right:
% spm_max creates peak height areas regardless of clusters - 
% possible to exclude a cluster if another voxel w/i 10 mm in a DIFFERENT
% cluster has a higher Z-score.  Also, if voxels in excursion set peak region
% are > radius mm apart, then it'll pick both edges of the exc. set as maxima
% and draw regions around those. 
% see tor_get_spheres3 for fixes.

% -----------------------------------------------
% get input arguments if empty
% -----------------------------------------------
if isempty(clusters), error('Must enter a cluster vector'),end
if  isempty(radius), radius = input('Enter radius of spheres in mm: '),end

if  isempty(SPM) |  isempty(VOL)
		disp(['Height threshold should be ' num2str(clusters(1).threshold)])
		[SPM,VOL,xX] = spm_getSPM;  
end

cl_out = [];

% -----------------------------------------------
% define 'excursion sets' - or sub-clusters based on local maxima
% -----------------------------------------------
[rsize rZ rXYZ rCl] = spm_max(SPM.Z,SPM.XYZ); 

% -----------------------------------------------
% exclude local maxima within 2 * radius of another local maximum
% -----------------------------------------------
rXYZ = cluster_sphere_limit_distance(rXYZ,2*radius,VOL.M,rZ);

% -----------------------------------------------
% define spheres around each local max
% uses find_sphere.m by Russ Poldrack, in spm roi toolbox
% needs mm coordinates.
% -----------------------------------------------
XYZmm = voxel2mm(rXYZ,VOL.M);

for i = 1:size(rXYZ,2)
	sphere(i)=find_sphere(XYZmm(:,i),radius,SPM,VOL);
end

% -----------------------------------------------
% find index of which sub-cluster sphere each voxel belongs to.
% -----------------------------------------------

for i = 1:length(clusters)
	clusters(i).subindex = sub_cluster_index(clusters(i),sphere);
end

% -----------------------------------------------
% get cluster information from existing clusters
% -----------------------------------------------

for i = 1:length(clusters)

	% get the cluster this sub-vol came from
	% -------------------------------------------
	cl = clusters(i);	

	myspheres = unique(cl.subindex);
	myspheres(myspheres == 0) = [];
	for j = myspheres
		cl_out(j).XYZmm = cl.XYZmm(:,cl.subindex == j);
		cl_out(j).XYZ = cl.XYZ(:,cl.subindex == j);
		cl_out(j).Z = cl.Z(cl.subindex == j);

		% get the timeseries, if timeseries data is saved.
		% -------------------------------------------
		if isfield(cl,'all_data')
			cl_out(j).all_data = cl.all_data(:,cl.subindex == j);
			cl_out(j).timeseries = nanmean(cl_out(j).all_data')';
		end
	
        cl_out(j).title = clusters(i).title;
        cl_out(j).threshold = clusters(i).threshold;
        cl_out(j).voxSize = clusters(i).voxSize;
        cl_out(j).title = clusters(i).title;
        
        if isfield(cl,'M'),cl_out(j).M = cl.M;,end
		cl_out(j).mm_center = center_of_mass(cl_out(j).XYZmm,cl_out(j).Z); %sphere(j).center';
        cl_out(j).center = mm2voxel(cl_out(j).mm_center,cl_out(j));
		cl_out(j).numVox = size(cl_out(j).XYZ,2);
		cl_out(j).radius = radius;
		cl_out(j).from_cluster = i;
		cl_out(j).name = [cl.name '_sph' num2str(j)];
        
        
	end
end

% ----------------------------------------------------------------------------------
% get the SNR for the cluster if it looks like rfx data rather than individual data
% ----------------------------------------------------------------------------------
for i = 1:length(cl_out)
    if length(cl_out(i).timeseries) < 60
        try
            cl_out(i).avg_snr = get_snr(cl_out(i).timeseries);
            cl_out(i).snr = get_snr(cl_out(i).all_data);
            fprintf(1,'sub-cluster %3.0f : avg SNR = %3.2f, range %3.2f - %3.2f \n',i,cl_out(i).avg_snr,min(cl_out(i).snr),max(cl_out(i).snr))
        catch
            disp(['tor_get_spheres2: problem getting SNR.'])
        end
    end
end
            
return
