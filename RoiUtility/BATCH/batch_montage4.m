function batch_montage4(varargin)
% run in the dir above results
% works on rfx, snpm, or correl directories
% whatever - looks for clusters files in subdirs.
% does surfaces as well.
%
% if varargin is entered, looks for t_corr_0002 clusters 
% rather than all clusters
% and changes output filenames slightly.
% input is (flag 1 or 0 to do this)
%

if length(varargin) > 0, 
    tflag = varargin{1};, 
    d = get_filename2('*/t_corr0002_clusters.mat');
else, 
    tflag = 0;
    d = get_filename2('*/*clusters.mat'); 
end

diary cluster_tables.txt

for i = 1:size(d,1), 
    
    [dd,ff] = fileparts(d(i,:));
    
    if tflag, ff = [dd '_t_corr_0002'];,end
    
    
    clear clusters; load(d(i,:));
    
    if tflag
        montage_clusters([],clusters,[2 2]),
        saveas(gcf,[ff '_montage_colbar'],'fig'), close
    else
        montage_clusters([],clusters),
    end

    diary on
    disp(ff)
    cluster_table(clusters);
    fprintf(1,'\n')
    diary off
    
    saveas(gcf,[ff '_montage'],'fig'),
    saveas(gcf,[ff '_montage'],'tif'), close
    
    dosurf = 0;
    
    if dosurf
        cluster_surf(clusters,which('surf_brain_render_T1.mat'),5,{[1 0 0]},'heatmap')
        set(gcf,'Color','w')
        saveas(gcf,[ff '_surf'],'fig')
        saveas(gcf,[ff '_surf'],'tif'), close
    
        saveas(gcf,[ff '_surf_colbar'],'tif'), close
    end

end

return
