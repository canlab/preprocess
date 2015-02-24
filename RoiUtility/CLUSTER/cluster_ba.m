function [clusters,strs,med2,um2,allcnts] = cluster_ba(clusters,whichclust)
% function [clusters,strs,clusvec,all_bas,ba_counts] = cluster_ba(clusters,whichclust)
%
% Takes a clusters structure and a cluster number, and finds the composition of voxels
% in that cluster (using the XYZmm field) according to Brodmann Areas and some labeled
% subcortical regions (thalamus, cerebellum, amygdala, hippocampus, caudate, putamen, etc.)
% Appends a string describing the composition to the clusters structure. 
%
% if whichclust is empty, loops through all clusters
% 
%
% by Tor Wager

if isempty(whichclust), whichclust = 1:length(clusters);,end
strs = [];

for cl = whichclust
    
    clear med2
    [medBA,tlab] = mni2BA(clusters(cl).XYZmm);

    for j=1:size(medBA,1);
        med2{j}=deblank(medBA(j,:));,
        if isempty(med2{j}) | strcmp(med2{j},'0'),
            med2{j} = deblank(tlab(j,:));,
        else
            med2{j} = ['BA' med2{j}];
        end
        if isempty(med2{j}),med2{j} = 'Unknown';,end
    end

        clear cnts
        um2 = unique(med2); for j=1:length(um2),cnts(j)=sum(strcmp(med2,um2{j}));,end
        allcnts = cnts;
        thr = max(round(sum(cnts)*.05),1);
        str = [];
        tot = sum(cnts);
        
        while any(cnts >= thr)
            
            wh = find(cnts == max(cnts));
            for i = 1:length(wh)
                
                s = sprintf('%s(%d%%), ',um2{wh(i)},round(100 * (cnts(wh(i))./tot)));
                str = [str s];
            end
            
            cnts(wh) = 0;
            
        end
        
        str = str(1:end-2);
        strs = str2mat(strs,str);
        clusters(cl).BAstring = str;
        
        
        
end
        
strs = strs(2:end,:);

return