function tmpcl2 = redefine_clusters(clpos)
% tmpcl2 = redefine_clusters(clpos)
% used in cluster_ttest.m, clusters_sigregions2.m
% tor wager

go=[]; for i=1:length(clpos),
    if isempty(clpos(i).XYZ),go(i)=0;,
    else,go(i)=1;, 
        if size(clpos(i).XYZmm,2) > 1,
            clpos(i).mm_center = mean(clpos(i).XYZmm');,
            clpos(i).numVox = size(clpos(i).XYZmm,2);
        else
            clpos(i).mm_center = clpos(i).XYZmm';
        end
    end,
end,
tmpcl=clpos(find(go));
tmpcl2 = [];
for i=1:length(tmpcl)
    for j = 1:max(tmpcl(i).clusters)
        tmp2(j) = tmpcl(i);
        tmp2(j).XYZ = tmpcl(i).XYZ(:,tmpcl(i).clusters == j);
        tmp2(j).XYZmm = tmpcl(i).XYZmm(:,tmpcl(i).clusters == j);
        tmp2(j).Z = tmpcl(i).Z(tmpcl(i).clusters == j);
        tmp2(j).t = tmpcl(i).t(:,tmpcl(i).clusters == j);
        tmp2(j).numVox = size(tmp2(j).XYZmm,2);
        tmp2(j).mm_center = mean(tmp2(j).XYZmm');,
        if length(tmp2(j).mm_center) == 1, tmp2(j).mm_center = tmp2(j).XYZmm';,end
    end
    tmpcl2 = [tmpcl2 tmp2];
end

return
