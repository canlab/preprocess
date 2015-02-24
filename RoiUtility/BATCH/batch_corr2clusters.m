function batch_corr2clusters(varargin)
% batch_corr2clusters
%
% by Tor Wager
%
% start this in the directory above each random effects/spm contrasts directory
% containing spm rfx analyses
%
% the script gets structures and clusters, and correlates each cluster
% with behavioral interference.

d = dir;
cwd = pwd;
diary([cwd filesep 'corr_cluster_summary.txt'])
diary off

if length(varargin) > 0, EXPT = varargin{1};
else
    
    if(exist('EXPT.mat') == 2), load EXPT, else, error('No EXPT.mat file.  Create with get_expt_info.m'),end
end

for i = 3:length(d)

    if d(i).isdir & strcmp(d(i).name(1:4),'corr'),eval(['cd ' d(i).name]),end
  
    if (exist('SPM.mat') == 2)
    
    
    % --------------------------------------------------------------------------------------------------
    % * locate contrast names and get contrast values
    % --------------------------------------------------------------------------------------------------
    
    EXPT.CORR.maskWOtherContrasts = 0; ...mask with other contrasts?
    EXPT.CORR.resdir = pwd;	...results directory

    % which image files (P) in EXPT.CORR? 
    %a = pwd; a = str2num(a(end-3:end));
    a = pwd; [dummy,a,dummy] = fileparts(a);
    switch(a(end-1:end))
    case 'ht', mypar = 1;
    case 'ay', mypar = 2;
    case 'pt', mypar = 3;
    otherwise, mypar = 0;, % disp(['directory name ' pwd ' is not height, delay, or intercept.'])
    end
    a = str2num(a(5:8));

    if mypar == 0
        % regular snpm directory
        P = EXPT.CORR.P{EXPT.CORR.connums == a};
    else
        EXPT.DX.connums = 1:length(EXPT.DX.connames);
        P = EXPT.DX.nlconP{mypar}{EXPT.DX.connums == a};
    end
    
    for myIc = 2:3    % loop thru contrasts - pos and neg
        
    EXPT.CORR.Ic = myIc;
    %disp(['Ic = ' num2str(EXPT.CORR.Ic)])
    [CLU,clusters] = spm2struct(pwd,EXPT.CORR,P);
    %[clusters] = tor_extract_rois(P,CLU,CLU,1);
    
    for j = 1:length(clusters)
        clusters(j).title = [clusters(j).title ' Ic = ' num2str(myIc)];
    end
    
    if length(clusters) > 100
        warning('More than 100 clusters!  Skipping save and table.')
        close all
    else
    % --------------------------------------------------------------------------------------------------
    % * correlate contrasts with behavioral vectors
    % --------------------------------------------------------------------------------------------------
    
    if isfield(EXPT,'behavior') 
        if ~isempty(EXPT.behavior) & isfield(clusters,'timeseries')
            clusters = simple_correl_clusters(clusters,EXPT.behavior,0);
        end
    else
        behav = input('Enter vector of numbers for behavioral scores, or return to skip correlations.');
        if ~isempty(behav)
            clusters = simple_correl_clusters(clusters,behav,1);
            EXPT.behavior = behav;
            save ../EXPT EXPT
        end
    end
    
    % --------------------------------------------------------------------------------------------------
    % * print table and save contrast clusters
    % --------------------------------------------------------------------------------------------------
    diary on
    fprintf(1,'\n')
    disp(['directory is ' pwd])
    if mypar == 0
       conname = deblank(EXPT.CORR.connames(EXPT.CORR.connums == a,:)); 
    else  
       conname = str2mat(EXPT.DX.connames);
       conname = deblank(conname(EXPT.DX.connums == a,:));
    end
    fprintf(1,'Contrast: %s',conname)
    cluster_table(clusters);
    diary off
    
    for j = 1:length(clusters), clusters(j).title = conname;,end
     a2 = pwd; [dummy,a2,dummy] = fileparts(a2);
    if myIc == 2
        eval(['save ' a2 '_pos_clusters clusters CLU'])   
    elseif myIc == 3
        eval(['save ' a2 '_neg_clusters clusters CLU']) 
    else
        error('Ic must be 2 or 3')
    end
    
    end % end more than 100 clusters
    
 end    % loop thru pos and neg contrasts
 
    % this way lets us choose t-maps for the correlations
    %[CLU,clusters] = snpm2struct(pwd,1);
    %a = pwd; a = a(end-7:end);
    %eval(['load ' a '_b4mask_trimrt_clusters'])
    
    % to get max fixed-effects z from *_b4data_b4mask_clusters
    % hack for including this column in the table for biman4
    %for j = 1:length(clusters)
    %    fprintf(1,'\n%3.2f',spm_t2z(max(mean(clusters(j).all_data) * sqrt(7)),7)  )
    %end
    %disp('  ')
    
   
    % get sub-clusters and correlations with sub-clusters
    %[sub_clusters,clusters] = tor_get_spheres2(clusters,10,CLU,CLU);
    %[sub_clusters,cmatx] = biman5_correl_clusters2(sub_clusters);   
    %a = pwd; a = a(end-7:end)
    % eval(['save ' a '_b4data_b4mask_tscore_clusters CLU clusters'])

    %O.which_cl = 1:length(clusters); O.clcol = []; O.dobrain = 'n';
    % renderCluster_ui(O)biman5
    
    end

    if d(i).isdir & strcmp(d(i).name(1:4),'corr'),cd .., end
    
end