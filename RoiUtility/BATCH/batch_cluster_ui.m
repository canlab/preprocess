% Tor Wager, 10/28/01

disp('This script extracts and analyzes all clusters for a single subject')
disp('and analyzes them with the design matrix in the SPM.mat file you select.')
disp('It also fits a deconvolution matrix to the data, if desired.')
disp('It uses the same SPM.mat file to define clusters and analyze the timeseries.')

% --------------------------------------------------------------------------
% * Get the files and SPM.mat
% --------------------------------------------------------------------------

disp('Select SPM analysis output with design matrix and SPM.mat')
[SPM,VOL,xX,xCon,xSDM] = spm_getSPM;

nsess = input('How many sessions?');

for i = 1:nsess
	fnames{i} = spm_get([1 1000],'*.img',['Select scans for session ' num2str(i)],pwd,0);
end

% --------------------------------------------------------------------------
% * Extract cluster data
% --------------------------------------------------------------------------

[clusters, SPM, xX, xCon] = tor_extract_rois(fnames,SPM,VOL,xX)


clusters = analyze_cluster_rois(clusters,xX,82);
disp('Cluster structure is in workspace as variable ''clusters''')


% --------------------------------------------------------------------------
% * Deconvolution
% --------------------------------------------------------------------------

% Inputs
% --------------------------------------------------------------------------

dodeconv = input('Do deconvolution? 1 or 0: ');
if dodeconv

	eres = length(Sess{1}.sf{1}) ./ (size(xX.X,1) ./ nsess);
	disp(['Number of samples in Sess.sf (time bins) per TR is calculated at ' num2str(eres) ' (SPM default is 16)'])

	tp = input('Number of time points (TRs) to deconvolve over? ');

	% Get deconvolution matrix
	% --------------------------------------------------------------------------

	[individualSessionsDX,shortSf] = tor_make_deconv_mtx2(evts_sf,tp,eres);
	totalSf = zeros(size(shorfSf{1},1),1);
	for i = 1:nsess
		totalSf = totalSf + shortSf{i};
	end
	[DX,totalSfout] = tor_make_deconv_mtx2(totalSf,tp,1);


	clusters = batch_deconv_clusters(clusters,individualSessionsDX,DX,tp);
 
	for i = 1:length(clusters)
    		figure; plot_cluster_hrfs2(clusters(i),nsess);
	end
end

a = [SPM.swd filesep SPM.title ' Clusters'];
b = SPM.title; b(b == ' ') = '_';
disp(['saving: ' b '.mat])
eval('save ' b ' clusters individualSessionsDX DX shorfSf totalSf')






