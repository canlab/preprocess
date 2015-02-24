function batch_snpm2clusters(varargin)
% batch_snpm2clusters
%
% by Tor Wager
%
% start this in the directory above each random effects/snpm contrasts directory
% containing snpm analyses
%
% the script gets structures and clusters, and correlates each cluster
% with behavioral interference.

mypwd = pwd;
diary([mypwd filesep 'cluster_summary.txt'])
diary off

d = dir('snpm*');

if length(varargin) > 0
    EXPT = varargin{1};
else
    
    if(exist(['.' filesep 'EXPT.mat']) == 2), load EXPT, 
    else, 
	    expfile = spm_get(1,'EXPT.mat','Choose EXPT.mat file');
	    eval(['load ' expfile])
    end
end

if isfield(EXPT,'behavior')
    doplot = 0; %input('Draw correlation scatterplots? (2/1/0), 2 is label, 1 is yes, 0 is no. ');
end

for i = 1:length(d)

    if d(i).isdir,eval(['cd(''' d(i).name ''')']), end
    
    if (exist('SnPMt_filtered.img') == 2) | (exist('SnPMt_neg_filtered.img') == 2)
    
    % get names of P from EXPT - only uses if the ones stored in the SnPM analyses don't work 
    % (e.g., if you move the files to a new computer.)
    
    myp = str2num(d(i).name(5:8));
    myp = myp == EXPT.SNPM.connums;
    myp = EXPT.SNPM.P{find(myp)};
    
    try
    [CLU,clusters] = snpm2struct(pwd,0,myp);
    catch
    end

    try
    [CLU_neg,clusters_neg] = snpm2struct_neg(pwd,0,myp);
    catch
    end
    
    % --------------------------------------------------------------------------------------------------
    % * locate contrast names and get contrast values
    % --------------------------------------------------------------------------------------------------
    
    %a = pwd; a = str2num(a(end-3:end));
    a = pwd; [dummy,a,dummy] = fileparts(a);
    switch(a(end-1:end))
    case 'ht', mypar = 1;
    case 'ay', mypar = 2;
    case 'pt', mypar = 3;
    otherwise, mypar = 0;, disp(['directory name ' pwd ' is not height, delay, or intercept.'])
    end
    a = str2num(a(5:8));
    
    if mypar == 0
        % regular snpm directory
        P = EXPT.SNPM.P{EXPT.SNPM.connums == a};

    else
        % nonlin directory
        EXPT.DX.connums = 1:length(EXPT.DX.connames);
        P = EXPT.DX.nlconP{mypar}{EXPT.DX.connums == a};      
    end
    
    if exist('CLU') == 1
        [clusters] = tor_extract_rois(P,CLU,CLU,1);
    end
    
    if exist('CLU_neg') == 1
        [clusters_neg] = tor_extract_rois(P,CLU,CLU,1);
    end
    %[subclusters] = tor_get_spheres3(clusters);
    
    % --------------------------------------------------------------------------------------------------
    % * correlate contrasts with behavioral vectors
    % --------------------------------------------------------------------------------------------------
    
    if isfield(EXPT,'behavior')
        if doplot == 2
            clusters = simple_correl_clusters(clusters,EXPT.behavior,EXPT.subjects);
            %subclusters = simple_correl_clusters(subclusters,EXPT.behavior,EXPT.subjects);
        elseif doplot == 1
            clusters = simple_correl_clusters(clusters,EXPT.behavior,1);
            %subclusters = simple_correl_clusters(subclusters,EXPT.behavior,1);
        else
                if exist('clusters') == 1, 
                    clusters = simple_correl_clusters(clusters,EXPT.behavior,0);
                end
            %subclusters = simple_correl_clusters(subclusters,EXPT.behavior,0);
        end
        
        if exist('clusters_neg') == 1
            clusters_neg = simple_correl_clusters(clusters_neg,EXPT.behavior,0);
        end
        
    else
        behav = input('Enter vector of numbers for behavioral scores, or return to skip correlations.');
        if ~isempty(behav)
            clusters = simple_correl_clusters(clusters,behav,EXPT.subjects);
            subclusters = simple_correl_clusters(clusters,behav,EXPT.subjects);
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
       conname = deblank(EXPT.SNPM.connames(EXPT.SNPM.connums == a,:)); 
    else  
       conname = str2mat(EXPT.DX.connames);
       conname = deblank(conname(EXPT.DX.connums == a,:));
    end

    fprintf(1,'Contrast: %s',conname)
    %cluster_table(clusters,subclusters)
    if exist('clusters') == 1
        disp('Positive effects')
        cluster_table(clusters);
            for j = 1:length(clusters), clusters(j).title = conname;,end
    end
    
    if exist('clusters_neg') == 1
        disp('Negative effects')
        cluster_table(clusters_neg);
            for j = 1:length(clusters_neg), clusters_neg(j).title = conname;,end
    end
    

     %a = pwd; a = a(end-7:end);
     a = pwd; [dummy,a,dummy] = fileparts(a);
     eval(['save ' a '_clusters clusters* CLU*'])   
    
     diary off
     
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

    if d(i).isdir,cd .., end
    
end