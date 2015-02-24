% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
%
%
% Routine to make a superclusters
%
% function [results, nClusters] = SOM_SuperCluster(SOM,nClusters)
%
% Grid must be square!
%
% This will return a series of supercluster suggestions
% The code will run through dot product thresholds from 
% 0.01 to 1.00
%
% Condition for SOM examplars combining into supercluster:
%
%    U.V >= Threshold
%    and nearest neighbors in the square grid.
%    
% The returned list of supercluster suggestions is based on 
% supercluster maps being different. It runs through the full
% list linearly, so there is a chance to get some discontinuity
% in the output - but most likely not!
%
% results   is a nGrid x nGrid x 100 matrix.
% nClusters is how many clusters are present in each supercluster
%              map.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [results, nClusters] = SOM_SuperCluster(SOM)

nGrid = sqrt(size(SOM,2));

% Make 100 super cluster maps!

results = zeros(nGrid,nGrid,100);

nClusters = [];

% Make an array of indexes. 

xs = (1:nGrid)'*ones(1,nGrid);
ys = ones(nGrid,1)*(1:nGrid);

% Make a second set of arrays for indexing the big cross array of SOM.

XS = (1:nGrid^2)'*ones(1,nGrid^2);
YS = ones(nGrid^2,1)*(1:nGrid^2);

% Call the routine that returns a pattern for nearest neighbors for
% any given element of the SOM.

nearestMatrix = SOM_Nearest(nGrid);

% Create a huge array of all U.V products between i,j pairs 
% in the SOM.

SOMbySOM = reshape(SOM'*SOM,nGrid*ones(1,4));

% Now make the generic super cluster map.

superCluster = reshape([1:nGrid*nGrid],[nGrid nGrid]);

% Find the index of the neighbors to look up the SOMbySOM value.

nearestIndex = find(nearestMatrix);

% Grab the U.V values of the closest neighbors.

UVValues = SOMbySOM(nearestIndex);

somData = [nearestIndex UVValues];

sortedData = sortrows(somData,2);

% Find out which SOM pairs are correlated.

SOMIdx1 = XS(sortedData(:,1));
SOMIdx2 = YS(sortedData(:,1));

% Get the index of each SOM on the the square grid.

% idx1X = xs(SOMIdx1);
% idx1Y = ys(SOMIdx1);
% 
% idx2X = xs(SOMIdx2);
% idx2Y = ys(SOMIdx2);

% diagMask = triu(1-diag(ones(1,nGrid^2)));
% 
% somValidx = find(diagMask);
% 
% CrossValues = SOMbySOM(somValidx);
% 
% Find the stats of a random map that should have no correlation in it.

% RandomMap = rand(size(SOM))-.5;
% 
% RandomMap = unitNormMatrix(RandomMap,1);
% 
% RByR = RandomMap'*RandomMap;
% 
% RRCV = RByR(somValidx);
% 

% Now collapse the supercluster maps.

for iThresh = 1:100
    pairThreshold = iThresh/100;
    cSuperCluster = superCluster;
    nMatches = 0 ;
    for iPair = size(sortedData,1):-1:1
        if sortedData(iPair,2) >= pairThreshold
            nMatches = nMatches + 1;
            findingIndex = 1;
            theIndex = SOMIdx1(iPair);
            while findingIndex==1
                if cSuperCluster(theIndex) == theIndex
                    findingIndex = 0;
                else    
                    theIndex = cSuperCluster(theIndex);
                end
                cSuperCluster(SOMIdx2(iPair)) = theIndex;
            end
        end
    end
    results(:,:,iThresh) = cSuperCluster;
    uniqueSOM = zeros(size(cSuperCluster));
    uniqueSOM(reshape(cSuperCluster,[1 prod(size(cSuperCluster))])) = 1;
    nClusters = [nClusters length(find(uniqueSOM))];
end

% Now pack the results.

superClusters = results(:,:,1);
nC = nClusters(1);
iClusterMaps = 1;
for iC = 2:100
    if any(any(results(:,:,iC)-results(:,:,iC-1)))
        nC = [nC nClusters(iC)];
        iClusterMaps = iClusterMaps + 1;
        superClusters(:,:,iClusterMaps) = results(:,:,iC);
    end
end

% Now renumber the SuperClusters in a linear fashion.

for iC = 1:size(superClusters,3)
    clusterNumbers = unique(superClusters(:,:,iC));
    sizeCluster = [];
    for iN = 1:length(clusterNumbers)
        sizeCluster = [sizeCluster prod(size(find(superClusters(:,:,iC)==clusterNumbers(iN))))];
    end            
    tmpArray = [clusterNumbers sizeCluster'];
    tmpArray = sortrows(tmpArray,2);
    clusterNumbers = flipud(tmpArray(:,1));
    sc = zeros(size(superClusters(:,:,iC)));
    for iN = 1:length(clusterNumbers)
        sc(find(superClusters(:,:,iC)==clusterNumbers(iN))) = iN;
    end
    superClusters(:,:,iC) = sc;
end

results = superClusters;
nClusters = nC;

%
% All done.
%
