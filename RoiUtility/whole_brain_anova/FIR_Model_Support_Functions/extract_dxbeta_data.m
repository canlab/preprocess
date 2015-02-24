function [clusters,EXPT] = extract_dxbeta_data(EXPT,clusters,varargin)
% [clusters,EXPT] = extract_dxbeta_data(EXPT,clusters,[plot])
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
% EXPT.FIR.dxbetas     % see group_fit_model and group_brain_nlfit
%                       % see also get_nlcon_names_for_rfx
% EXPT.FIR.dxtrialonsets % see tor_build_deconv_design_ui
% EXPT.FIR.numframes
% EXPT.FIR.dxnames       % cell array of trial type names
%
% For plotting:
% EXPT.FIR.TR        % repetition time of acquisition
%
% Last modified 3/3/05 by tor wager
%
% see also GROUP_BRAIN_NLFIT.M get_nlcon_names_for_rfx tor_build_deconv_design_ui
% check_cluster_data.m plot_dx_hrfs.m nl_hrf_plot.m
%
% Other stuff you might want to define:
% EXPT.FIR.mcol = {'ro-' 'go-' 'bo-' 'co-' 'mo-'}
% EXPT.FIR.indiv = 1;
% This requires EXPT.beh --> DOES MEDIAN SPLIT ON FIRST COLUMN
% EXPT.FIR.regsofinterest = 1:5;
% EXPT.FIR.baseline = 1:2;
%
% After extracting and saving in HRF substructure, to plot:
%     EXPT = plot_dx_hrfs(EXPT,clusters);
%     EXPT = plot_dx_hrfs_indiffs(EXPT,clusters);

% DIAGNOSIS
% cl{1}(5).all_data(1:5,1:5)
%EXPT.subjects
%cl{1}(5).XYZ(:,1:5),cl{1}(5).XYZmm(:,1:5)
%

%doplot = input('Plot HRFs (1/0)?');
doplot = 0;
if length(varargin) > 0, doplot = varargin{1};,end


EXPT = CreateExpt('extractdx');

for i = 1:length(EXPT.subjects)
    
    disp('  ')
    disp(['Starting subject ' num2str(i)])
    
    % extract each cluster separately to preserve orig clusters
    %for j = 1:length(clusters)
        
        %CLU = clusters2CLU(clusters(j),clusters(1).M);
        
        % WARN: if vox sizes or origin is not the same, must transform here if different
        V = spm_vol(deblank(EXPT.FIR.dxbetas{i}(1,:)));
        if any(any(clusters(1).M(1:3,1:4) - V.mat(1:3,1:4)))
            for j = 1:length(clusters)
                CLU = clusters2CLU(clusters(j),clusters(1).M);
                fprintf(1,'\t(Transforming cluster %3.0f voxel coords to beta img space)\n',j)
                CLU = transform_coordinates(CLU,V.mat);
                clusters(j).M = V.mat;clusters(j).VOX = diag(V.mat(1:3,1:3))';
                clusters(j).voxSize = diag(V.mat(1:3,1:3))';
                clusters(j).XYZ = CLU.XYZ; 
                clusters(j).XYZmm = CLU.XYZmm;
                clusters(j).Z = ones(1,size(clusters(j).XYZ,2));
            end
        end
          
        cl{i} = tor_extract_rois(EXPT.FIR.dxbetas{i},clusters);
        
        for j = 1:length(clusters)
        
            % index so 1st is cluster, matrix of subjects
            ts{j}(:,i) = cl{i}(j).timeseries;
        
        end

end

disp(['Getting group HRFs'])

st = cumsum([1 EXPT.FIR.numframes]);   
en = st(2:end) - 1;         % ending values
st = st(1:end-1);           % starting values


for j = 1:length(clusters)
    
    clusters(j).HRF.indiv = ts{j};
    clusters(j).HRF.avg = nanmean(ts{j}')';
    
    % now break into HRFs
    %for i = 1:length(EXPT.FIR.dxtrialonsets-1)
    %    tmp = clusters(j).HRF.avg(EXPT.FIR.dxtrialonsets(i):EXPT.FIR.dxtrialonsets(i)+EXPT.FIR.numframes-1);
    
    for i = 1:length(st) % for each condition -- new way
        tmp = clusters(j).HRF.avg(st(i):en(i));    
        
        eval(['clusters(j).HRF.' EXPT.FIR.dxnames{i} ' = tmp'';'])
    end

    %clusters(j).HRF = HRF;

end

for j = 1:length(clusters)
    for i = 1:size(clusters(j).HRF.indiv,2)
        
        %m = beta2conditions(clusters(j).HRF.indiv(:,i)',EXPT.FIR);
        for k = 1:length(st) % for each condition -- new way
            m{k} = clusters(j).HRF.indiv(st(k):en(k),i);        % get this ind hrf for this subj 
        end
        
        
        for k = 1:length(m)
            b(i,k) = mean(m{k}(EXPT.FIR.baseline));  % estimate mean baseline
            m{k} = m{k} - b(i,k);                   % and subtract it
            
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
        


% individual diffs

if EXPT.FIR.indiv
    
    try, beh = EXPT.beh;  , 
    
    catch, 
        try
            beh = EXPT.cov;
        catch
            beh = EXPT.behavior';, 
        end
    end
    
    % take only first column, if multiple columns
    if size(beh,2) == 2, beh = beh(:,1);,end
    
low = find(beh < median(beh));
high = find(beh > median(beh));
    
for j = 1:length(clusters)
    
    lowdat = clusters(j).HRF.indiv(:,low);
    hidat = clusters(j).HRF.indiv(:,high);
    
    
    % lows
    for i = 1:size(lowdat,2)
        
        for k = 1:length(st) % for each condition -- new way
            m{k} = lowdat(st(k):en(k),i);        % get this ind hrf for this subj 

            b(i,k) = mean(m{k}(EXPT.FIR.baseline));  % estimate mean baseline
            m{k} = m{k} - b(i,k);                   % and subtract it
            
            clusters(j).HRF.lowhrf{k}(i,:) = m{k};
            clusters(j).HRF.lowbase(i,k) = b(i,k);
        end
    end
    
    for k = 1:length(m)
            clusters(j).HRF.lowHRF{k} = mean(clusters(j).HRF.lowhrf{k});
            clusters(j).HRF.lowSTE{k} = ste(clusters(j).HRF.lowhrf{k});
    end
        
    clusters(j).HRF.lowBASE = mean(clusters(j).HRF.lowbase);
    clusters(j).HRF.lowBASESTE = ste(clusters(j).HRF.lowbase);
    
    
    % highs
    for i = 1:size(hidat,2)
        
        for k = 1:length(st) % for each condition -- new way
            m{k} = hidat(st(k):en(k),i);        % get this ind hrf for this subj 

            b(i,k) = mean(m{k}(EXPT.FIR.baseline));  % estimate mean baseline
            m{k} = m{k} - b(i,k);                   % and subtract it
            
            clusters(j).HRF.hihrf{k}(i,:) = m{k};
            clusters(j).HRF.hibase(i,k) = b(i,k);
        end
    end
    
    for k = 1:length(m)
            clusters(j).HRF.hiHRF{k} = mean(clusters(j).HRF.hihrf{k});
            clusters(j).HRF.hiSTE{k} = ste(clusters(j).HRF.hihrf{k});
    end
        
    clusters(j).HRF.hiBASE = mean(clusters(j).HRF.hibase);
    clusters(j).HRF.hiBASESTE = ste(clusters(j).HRF.hibase);
    
end % clusters

end % if indiv



% plot option

if doplot
    if ~isfield(EXPT.FIR,'indiv'),EXPT.FIR.indiv = 0;, end
    
    EXPT = plot_dx_hrfs(EXPT,clusters,0,1,3,EXPT.FIR.indiv);
    
    %if indiv
        
    %    EXPT = plot_dx_hrfs_indiffs(EXPT,clusters,0,1,3);
        
    %end
end





return


