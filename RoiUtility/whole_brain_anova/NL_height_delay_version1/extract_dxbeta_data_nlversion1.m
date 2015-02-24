function clusters = extract_dxbeta_data(EXPT,clusters,varargin)
% clusters = extract_dxbeta_data(EXPT,clusters,[plot])
%
% tor wager
% Extracts HRF esatimate data for each subject in each cluster
% averaging over cluster voxels for each individual (saved in clusters.all_data)
% and averaging over subjects (saved in clusters.timeseries)
% 
% inputs:  EXPT, as defined with get_expt_info.m
%          clusters, as defined with tor_extract_rois.m
%          optional: plot flag for HRF averages 
%
% Requires fields:
%
% clusters.M        % SPM mat file for volume mapping
% other cluster fields defined in tor_extract_rois.m
%
% EXPT structure: see get_expt_info.m
% EXPT = tor_build_deconv_design_ui(EXPT);
%
% EXPT.dx_beta_imgs     % see group_fit_model and group_brain_nlfit
%                       % see also get_nlcon_names_for_rfx
% EXPT.DX.dxtrialonsets % see tor_build_deconv_design_ui
% EXPT.DX.numframes
% EXPT.DX.dxnames       % cell array of trial type names
%
% For plotting:
% EXPT.DX.TR        % repetition time of acquisition
%
% Last modified 3/20/04 by tor wager
%
% see also GROUP_BRAIN_NLFIT.M get_nlcon_names_for_rfx tor_build_deconv_design_ui
% check_cluster_data.m plot_dx_hrfs.m nl_hrf_plot.m

% DIAGNOSIS
% cl{1}(5).all_data(1:5,1:5)
%EXPT.subjects
%cl{1}(5).XYZ(:,1:5),cl{1}(5).XYZmm(:,1:5)
%

%doplot = input('Plot HRFs (1/0)?');
doplot = 0;
if length(varargin) > 0, doplot = varargin{1};,end

if doplot
    if ~isfield(EXPT.DX,'TR'), EXPT.DX.TR = input('Please enter TR in s: ');,end
end

for i = 1:length(EXPT.subjects)
    
    disp('  ')
    disp(['Starting subject ' num2str(i)])
    
    % extract each cluster separately to preserve orig clusters
    %for j = 1:length(clusters)
        
        %CLU = clusters2CLU(clusters(j),clusters(1).M);
        
        % WARN: if vox sizes or origin is not the same, must transform here if different
        V = spm_vol(deblank(EXPT.dx_beta_imgs{i}(1,:)));
        if any(any(clusters(1).M(1:3,1:4) - V.mat(1:3,1:4)))
            for j = 1:length(clusters)
                CLU = clusters2CLU(clusters(j),clusters(1).M);
                fprintf(1,'\t(Transforming cluster %3.0f voxel coords to beta img space)\n',j)
                CLU = transform_coordinates(CLU,V.mat);
                clusters(j).M = V.mat;clusters(j).VOX = diag(V.mat(1:3,1:3))';clusters(j).voxSize = diag(V.mat(1:3,1:3))';
                clusters(j).XYZ = CLU.XYZ; clusters(j).XYZmm = CLU.XYZmm;
            end
        end
          
        cl{i} = tor_extract_rois(EXPT.dx_beta_imgs{i},clusters);
        
        for j = 1:length(clusters)
        
            % index so 1st is cluster, matrix of subjects
            ts{j}(:,i) = cl{i}(j).timeseries;
        
        end

    end

end

disp(['Getting group HRFs'])
for j = 1:length(clusters)
    
    clusters(j).HRF.indiv = ts{j};
    clusters(j).HRF.avg = nanmean(ts{j}')';
    
    % now break into HRFs
    for i = 1:length(EXPT.DX.dxtrialonsets-1)
        tmp = clusters(j).HRF.avg(EXPT.DX.dxtrialonsets(i):EXPT.DX.dxtrialonsets(i)+EXPT.DX.numframes-1);
        eval(['clusters(j).HRF.' EXPT.DX.dxnames{i} ' = tmp'';'])
    end

    %clusters(j).HRF = HRF;

end

for j = 1:length(clusters)
    for i = 1:size(clusters(j).HRF.indiv,2)
        m = beta2conditions(clusters(j).HRF.indiv(:,i)',EXPT.DX);
        for k = 1:length(m)
            b(i,k) = mean(m{k}(1:2));
            m{k} = m{k} - b(i,k);
            
            clusters(j).HRF.hrf{k}(i,:) = m{k};
            clusters(j).HRF.base(i,k) = b(i,k);
        end
    end
    
    for k = 1:length(m)
            clusters(j).HRF.HRF{k} = mean(clusters(j).HRF.hrf{k});
            clusters(j).HRF.STE{k} = ste(clusters(j).HRF.hrf{k});
    end
        
    clusters(j).HRF.BASE = mean(clusters(j).HRF.base);
    clusters(j).HRF.BASESTE = ste(clusters(j).HRF.base);
end
        
% plot option

keyboard

if doplot
    EXPT = plot_dx_hrfs(EXPT,clusters);
    
end


return


