function [s,ss] = silhouetteplot(X,varargin)
% function [s,ss] = silhouetteplot(X,varargin)
%
% plots average and complete linkage silhouette values by number of
% clusters in solution; hierarchical clustering using clusterdata
%
% varargin suppresses plot creation

if size(X,1) < 3, disp('NO Sil plot possible: Not enough data.'), return, end

if length(varargin) == 0, [h,col] = tor_fig;, end

    for j = 1:min(15,size(X,1))
        T = clusterdata(X,'maxclust',j,'linkage','average');
        [s] = silhouette(X,T);
        S(j) = mean(s);
        V(j) = ste(s);
        
        T = clusterdata(X,'maxclust',j,'linkage','complete');
        [ss] = silhouette(X,T);
        SS(j) = mean(ss);
        VV(j) = ste(s);
        
        T = clusterdata(X,'maxclust',j,'linkage','complete');
        [sss] = silhouette(X,T);
        SSS(j) = mean(ss);
        VVV(j) = ste(s);
        
    end

    plot(S,'ks-','MarkerFaceColor',[.5 .5 .5]);
    plot(SS,'ko--','MarkerFaceColor',[.1 .1 .1]);
    plot(SSS,'kv:','MarkerFaceColor',[0 0 0]);
    tor_bar_steplot(S,V,{'k'}); tor_bar_steplot(SS,-VV,{'k'});
    %tor_bar_steplot(S,-V,{'k'}); tor_bar_steplot(SS,-VV,{'k'});
    
    set(gca,'XTick',1:16); xlabel('Number of clusters'), ylabel('Avg silhouette value')
    legend({'Average linkage' 'Complete linkage' 'Single linkage'})
    
    
    return
    