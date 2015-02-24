function batch_cluster_tables(u,varargin)
% batch_cluster_tables(u,[k],[EXPT])
%
% u is height threshold p-value
% k is cluster extent threshold
%
% tor wager

if length(varargin) > 1, EXPT = varargin{2};

else 
	if exist('EXPT.mat') == 2, load EXPT, else, error('No EXPT.mat file in current directory'), end
end

if length(varargin) > 0, k = varargin{1};, else, k = 0;, end

% ----------------------------------------------------------------------------------
% * rfx directories
% ----------------------------------------------------------------------------------

% remove old clusters.mat files to avoid confusion
disp('Removing old t_rfx0002_clusters.mat files')
d = dir('rfx00*'); d = d(cat(1,d.isdir));,d = str2mat(d.name);
for i = 1:size(d,1)    
    try, eval(['!del ' deblank(d(i,:)) filesep 't_rfx0002_clusters.mat']),catch,end
    try, eval(['!rm ' deblank(d(i,:)) filesep 't_rfx0002_clusters.mat']),catch,end
end

disp('thresholding SPM t images and extracting clusters')
df = size(EXPT.SNPM.P{1},1) - 1;
threshold_spm_t(u,df,k,'pos',EXPT.SNPM.P);

diary batch_cluster_tables_output.txt
diary off

d = dir('rfx00*'); d = d(cat(1,d.isdir));,d = str2mat(d.name);

for i = 1:size(d,1)
    
    try
        load([deblank(d(i,:)) filesep 't_rfx0002_clusters'])
        clusters = simple_correl_clusters(clusters,EXPT.behavior,0);
        eval(['save ' deblank(d(i,:)) filesep 't_rfx0002_clusters clusters'])
        diary on
        disp(' ')
        disp([deblank(d(i,:)) ' ' EXPT.SNPM.connames(i,:)])
        cluster_table(clusters);
        diary off
    catch
        diary on
        disp(' ')
        disp([deblank(d(i,:)) ' ' EXPT.SNPM.connames(i,:)])
        disp('No significant results.')
        diary off
    end
    
end

% ----------------------------------------------------------------------------------
% * corr directories
% ----------------------------------------------------------------------------------

% remove old clusters.mat files to avoid confusion
disp('Removing old t_corr0002_clusters.mat files')
d = dir('corr00*'); d = d(cat(1,d.isdir));,d = str2mat(d.name);
for i = 1:size(d,1)    
    try, eval(['!del ' deblank(d(i,:)) filesep 't_corr0002_clusters.mat']),catch,end
    try, eval(['!rm ' deblank(d(i,:)) filesep 't_corr0002_clusters.mat']),catch,end
end


disp('thresholding SPM t images for correlations and extracting clusters')
df = size(EXPT.CORR.P{1},1) - 2;
threshold_corr_t(u,df,k,'pos',EXPT.CORR.P);

d = dir('corr00*'); d = d(cat(1,d.isdir));,d = str2mat(d.name);

for i = 1:size(d,1)
    
    try
        load([deblank(d(i,:)) filesep 't_corr0002_clusters'])
        clusters = simple_correl_clusters(clusters,EXPT.behavior,0);
        eval(['save ' deblank(d(i,:)) filesep 't_corr0002_clusters clusters'])
        diary on
        disp(' ')
        disp([deblank(d(i,:)) ' ' EXPT.CORR.connames(i,:)])
        cluster_table(clusters);
        diary off
    catch
        diary on
        disp(' ')
        disp([deblank(d(i,:)) ' ' EXPT.CORR.connames(i,:)])
        disp('No significant results.')
        diary off
    end
    
end


% ----------------------------------------------------------------------------------
% * robust directories
% ----------------------------------------------------------------------------------

% remove old clusters.mat files to avoid confusion
disp('Removing old rob_*_clusters.mat files')
d = dir('robust00*'); d = d(cat(1,d.isdir));,d = str2mat(d.name);
for i = 1:size(d,1)    
    try, eval(['!del ' deblank(d(i,:)) filesep 'rob_*_clusters.mat']),catch,end
    try, eval(['!rm ' deblank(d(i,:)) filesep 'rob_*_clusters.mat']),catch,end
end


disp('thresholding ROBUST t images for correlations and extracting clusters')
df = size(EXPT.SNPM.P{1},1) - size(EXPT.cov,2) - 1;
threshold_robust_t(u,df,k,'pos',EXPT.SNPM.P);

d = dir('robust00*'); d = d(cat(1,d.isdir));,d = str2mat(d.name);

for i = 1:size(d,1)
    
	dd = dir([d(i).name filesep 'rob_tmap*img']);

	for j = 1:length(dd)

    try
        load([deblank(d(i,:)) filesep 'rob_' dd(j).name(10:13) '_clusters'])
        clusters = simple_correl_clusters(clusters,EXPT.cov(:,1),0);
        eval(['save ' deblank(d(i,:)) filesep 'rob_' dd(j).name(10:13) '_clusters clusters'])
        diary on
        disp(' ')
        disp([deblank(d(i,:)) ' ' EXPT.SNPM.connames(i,:)])
        cluster_table(clusters);
        diary off
    catch
        diary on
        disp(' ')
        disp([deblank(d(i,:)) ' ' EXPT.SNPM.connames(i,:)])
        disp('No significant results.')
        diary off
    end
    
	end	% loop through results within dir

end
return
    
    