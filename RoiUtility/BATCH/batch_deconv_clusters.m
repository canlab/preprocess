function clusters = batch_deconv_clusters(clusters,DX,oDX,tp)
% function batch_deconv_clusters(clusters,DX,oDX,tp)
%
% DX is deconvolution matrix with individual estimates for each session
% oDX is deconvolution irrespective of session

% adjustedy is computed with analyze_cluster_rois

fprintf(1,'cluster ')

for i = 1:length(clusters)
    
    fprintf(1,'%d ',i)
    
    clusters(i).Db = pinv(DX) * clusters(i).adjustedy;
    clusters(i).oDb = pinv(oDX) * clusters(i).adjustedy;
    
    index = 1;
    for j = 1:tp:length(clusters(i).oDb)
        clusters(i).hrf_est{index} = clusters(i).oDb(j:j+tp-1);

        index = index + 1;
    end
        
end

return