function clusters = define_cluster(SPM,VOL)

	pThresh = .1;		% for one-tailed test

    % ----------------------------------------------------------------------------------
	% get cluster index number from each voxel
	% ----------------------------------------------------------------------------------
	cl_index = spm_clusters(SPM.XYZ);
	
	% ----------------------------------------------------------------------------------
	% define each cluster as cell in array.
	% ----------------------------------------------------------------------------------
	clusters = [];
	for i = 1:max(cl_index)
		disp(['Extracting cluster ' num2str(i) ' of ' num2str(max(cl_index))])
		a = find(cl_index == i);

		cl.title = SPM.title;
		cl.threshold = SPM.u;
		cl.voxSize = VOL.VOX;
		cl.name = [cl.title '_' num2str(i) '_' mat2str(size(a,2)) '_voxels'];
		cl.numVox = size(a,2);
		cl.Z = SPM.Z(a);
		cl.XYZmm = SPM.XYZmm(:,a);
		cl.XYZ = SPM.XYZ(:,a);
		try
			cl.pVoxelLev = spm_P(1,0,max(cl.Z),SPM.df,SPM.STAT,VOL.R,SPM.n);
			cl.pClustLev = spm_P(1,cl.numVox/prod(VOL.FWHM),SPM.u,SPM.df,SPM.STAT,VOL.R,SPM.n);
		catch
			warning('Can''t get SPM p voxel and cluster values.  Skipping spm_P.')
		end
        
	if isfield(SPM,'t') & isfield(SPM,'p')
        	cl.t = SPM.t(a);
        	cl.p = SPM.p(a);
        	cl.con1_XYZ = cl.XYZ(:,cl.p <= pThresh & (cl.t > 0));
        	cl.con2_XYZ = cl.XYZ(:,cl.p <= pThresh & (cl.t < 0));
        	cl.fuzzy = cl.XYZ(:,cl.p > pThresh);

		[h,cl.cl_p,ci,cl.cl_t] = t_test2(cl.t(~isnan(cl.t)));
        end

	if isfield(SPM,'imnames')
		cl.imnames = SPM.imnames;
	end
 
	clusters = [clusters,cl];

    end
    
 return