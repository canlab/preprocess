function [cl_out] = tor_get_spheres3(clusters)
%  [cl_out] = tor_get_spheres3(clusters)
%
% function: to create a vector of clusters (cl_out) consisting of thresholded
%	activation clusters broken into sub-clusters based on local maxima.
%	Sub-cluster centers are defined with spm_max, and
%	cluster boundaries are defined by a sphere of a chosen radius, in mm.
%	
%	To be included in a sub-cluster, a voxel must be within the sphere
%	and must be significant at the chosen threshold.
%
% OLD COMMENTS...
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
%

% -----------------------------------------------
% get input arguments if empty
% -----------------------------------------------
if isempty(clusters), error('Must enter a cluster vector'),end
%if  isempty(radius), radius = input('Enter radius of spheres in mm: '),end

%if  isempty(SPM) |  isempty(VOL)
%		disp(['Height threshold should be ' num2str(clusters(1).threshold)])
%		[SPM,VOL,xX] = spm_getSPM;  
%end

cl_out = [];
cli = 1;

for i = 1:length(clusters)
    
   	% get the cluster this sub-vol came from
	% -------------------------------------------
	cl = clusters(i);
    
    % -----------------------------------------------
    % define 'excursion sets' - or sub-clusters based on local maxima
    % -----------------------------------------------
    [rsize rZ rXYZ rCl] = spm_max(cl.Z,cl.XYZ); 

    for j = 1:max(rCl)
        
        cl_out(cli).from_cluster = i;
        cl_out(cli).M = cl.M;
        cl_out(cli).voxSize = cl.voxSize;
        cl_out(cli).title = [cl.title ' subclust'];
        cl_out(cli).threshold = cl.threshold;
        
        cl_out(cli).numVox = length(rZ);
        cl_out(cli).XYZ = rXYZ(:,rCl == j);
        cl_out(cli).Z = rZ(rCl == j);
        cl_out(cli).XYZmm = voxel2mm(cl_out(cli).XYZ,cl_out(cli).M);
  		cl_out(cli).mm_center = center_of_mass(cl_out(cli).XYZmm,cl_out(cli).Z); %sphere(cli).center';
        cl_out(cli).center = mm2voxel(cl_out(cli).mm_center,cl_out(cli));      
        
        wvox = match_rows(cl_out(cli).XYZ',cl.XYZ');

        % get the timeseries, if timeseries data is saved.
		% -------------------------------------------
		if isfield(cl,'all_data')
			cl_out(cli).all_data = cl.all_data(:,wvox);
            if size(cl_out(cli).all_data,2) > 1
			    cl_out(cli).timeseries = nanmean(cl_out(cli).all_data')';
            else
                cl_out(cli).timeseries = cl_out(cli).all_data;
            end
        end
	
		cl_out(cli).name = [cl.name '_subc' num2str(j)];
        
        cli = cli + 1;
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
            %fprintf(1,'sub-cluster %3.0f : avg SNR = %3.2f, range %3.2f - %3.2f \n',i,cl_out(i).avg_snr,min(cl_out(i).snr),max(cl_out(i).snr))
        catch
            disp(['tor_get_spheres2: problem getting SNR.'])
        end
    end
end
            
return
