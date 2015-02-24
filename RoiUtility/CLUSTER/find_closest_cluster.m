function [clusters,wh] = find_closest_cluster(allcl,i)
    % function [clusters,wh] = find_closest_cluster(allcl,i)
    %
    % TWO MODES:
    % Mode 1:
    % input: allcl = a cell array where each cell is a clusters structure with multiple clusters
    %        i = the index of the cluster to match to, in the first cell of clusters
    %
    % output: a clusters structure containing
    %           - the indexed (ith) cluster of cell 1
    %           - the closest clusters from each other cell in allcl
    %         wh: vector of indices for which cluster is chosen
    %
    % Mode 2:
    % input:  allcl = a clusters structure
    %          i = a set of x, y, z mm coordinates (row vector)
    %
    % % see also overlap.m and overlap_closest.m



    % Mode 1
    if iscell(allcl)

        mm = mean(allcl{1}(i).XYZmm,2);
        clusters(1) = allcl{1}(i);
        wh(1) = i;

        for j = 2:length(allcl)

            clear d
            for k = 1:length(allcl{j})

                ctr = mean(allcl{j}(k).XYZmm,2);
                d(k) = sum((mm - ctr) .^2);

            end

            whj = find(d == min(d));
            wh(j) = whj(1);
            clusters(j) = allcl{j}(wh(j));

        end


    else
        % Mode 2

        % i is x y z MM coords

        if iscol(i), i = i'; end
        
        for c = 1:length(allcl)

            xyz = allcl(c).XYZmm(1:3, :)';  % columns

            d = repmat(i, size(xyz,1), 1) - xyz; % diffs

            d = sum(d .^ 2, 2);  % sum of sq distances for each voxel          

            sumsq(c) = min(d);  % ssq for closest voxel

        end

        wh = find(sumsq == min(sumsq));

        clusters = allcl(wh);



    end


    return

