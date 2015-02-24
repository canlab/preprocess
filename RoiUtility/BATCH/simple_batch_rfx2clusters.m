% simple batch rfx2clusters

d = dir('rfx*'); diary([pwd filesep 'rfx_cluster_summary005.txt']);
for i = 1:length(d)
    if d(i).isdir
        cd(d(i).name)
        dd=dir('rfx*clusters.mat');
        try
            load(dd.name)
            clusters = cluster_table(clusters);
        catch
        end
        cd ..
    end
end
diary off

d = dir('corr*'); diary([pwd filesep 'corr_cluster_summary005.txt']);
for i = 1:length(d)
    if d(i).isdir
        cd(d(i).name)
        dd=dir('corr*pos_clusters.mat');
        try
            load(dd.name)
            clusters = cluster_table(clusters);
        catch
        end
        cd ..
    end
end
diary off


d = dir('rfx*'); diary([pwd filesep 'rfx_cluster_summary001.txt']);
for i = 1:length(d)
    if d(i).isdir
        cd(d(i).name)
        dd=dir('rfx*clusters.mat');
        try
            load(dd.name)
            clusters = cluster_table(clusters);
        catch
        end
        cd ..
    end
end
diary off

d = dir('corr*'); diary([pwd filesep 'corr_cluster_summary001.txt']);
for i = 1:length(d)
    if d(i).isdir
        cd(d(i).name)
        dd=dir('corr*pos_clusters.mat');
        try
            load(dd.name)
            clusters = cluster_table(clusters);
        catch
        end
        cd ..
    end
end
diary off