% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
%
% function results = SOM_CalculateMap(theData,nSOM,nIter)
%
% theData = theData(nSpace,nTime);
% nSOM    = number of basis functions in the SOM
% nIter   = number of iterations.
%
% results = results structure, it will at least contain
%           the resulting SOM!.
%
% results = results.SOM  (the map)
%           results.IDX  (for each space element which SOM it is related to
%                         best)
%           results.WT   (how much of the variance is explained by that)
%
%   this is a work in progress!
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function results = SOM_CalculateMap(theData,nSOM,nIter)

global SOM

% Which method to determine neighborhood (same code, just cleaned
% up?).

if ~isfield(SOM,'OldMethod') 
  SOM.OldMethod = 0
end

% Check to see if the necessary fields exist.

% The neighborhood size
if ~isfield(SOM,'sigma')
    SOM.sigma = 4;
    fprintf('SOM.sigma                -> %f\n',SOM.sigma);
end

% How quickly to modify the neighborhood size.
if ~isfield(SOM,'sigmaTimeConstant')
    SOM.sigmaTimeConstant = 1/4;
    fprintf('SOM.sigmaTimeConstant    -> %f\n',SOM.sigmaTimeConstant);
end

% How quickly to modify the map learning rate.
if ~isfield(SOM,'learningTimeConstant')
    SOM.learningTimeConstant = 2;
    fprintf('SOM.learningTimeConstant -> %f\n',SOM.learningTimeConstant);
end

% The initial map learning rate.
if ~isfield(SOM,'alpha')
    SOM.alpha = .1;
    fprintf('SOM.alpha                -> %f\n',SOM.alpha);
end

% Determine the size of the space we are dealing with.

if sqrt(nSOM)~=floor(sqrt(nSOM))
    fprintf('Can only do square SOM!\n');
    fprintf('Forcing SOM to be square\n');
    nSOM = floor(sqrt(nSOM)+.5)^2;
end

nGrid = sqrt(nSOM);

nSpace = size(theData,1);
nTime  = size(theData,2);

% Randomize our initial SOM and normalize the basis.
SelfOMap = rand(nTime,nSOM)-.5;
SelfOMap = unitNormMatrix(SelfOMap,1);
clear unitNormMatrix;

% Unit norm the data. We are not necessarily interested in effect
% size at the moment versus explanatory basis. We can 
% get the "beta's" later!

theDataN = unitNormMatrix(theData,2);
clear unitNormMatrix;

% Ok, let's go iterate the map. At the moment
% we'll just stop the calculation after a sufficient number
% iterations.

iter = 0;

neighborDist = SOM_NeighborDist(nGrid);

lneighborDist = reshape(neighborDist,[nGrid nGrid nGrid*nGrid]);

while iter < nIter
    tic;
    % Increment the iteration
    iter = iter + 1;
    % Learning function = g(time);
    alpha = SOM_Alpha(iter,nIter);
    % Determine the neighbor map.
    % For each element of the SOM, there will be an associated
    % neighborhood. The call should return the indices of these
    % SOM vectors as well as the distance for use in the weight 
    % calculation.    
    NeighSigma = SOM.sigma*exp(-iter/nIter/SOM.sigmaTimeConstant);
    if SOM.OldMethod == 0
      SOMNeighborMap = exp(-lneighborDist.^2/NeighSigma^2);    
    else
      SOMNeighborMap = SOM_NeighborMap([nGrid nGrid],NeighSigma);    
    end
    % Determine the SOM vectors the data are closest to.
    idx = SOM_FindClosest(theDataN,SelfOMap);
    clear SOM_FindClosest;
    % Now calculate the perturbation to the SOM.
    dSelfOMap = 0*SelfOMap;
    nullSOM   = 0*SelfOMap;
    for iSOM = 1:nSOM
        % Deterimine which data affect this SOM vector.
        dI = find(idx==iSOM);
        if length(dI) > 0
            pV = sum(theDataN(dI,:),1)';
            if SOM.OldMethod == 0
              dSelfOMap = dSelfOMap + SOM_ModBasis(pV,nullSOM, ...
                                                   SOMNeighborMap(:,:,iSOM),iter,nIter);
            else
              dSelfOMap = dSelfOMap + SOM_ModBasis(pV,nullSOM, ...
                                                   SOMNeighborMap{iSOM},iter,nIter);
            end            
            clear SOM_ModBasis;
        end
    end
    % Update the SOM.
    SelfOMap = SelfOMap + alpha*dSelfOMap;
    SelfOMap = unitNormMatrix(SelfOMap,1);
    clear unitNormMatrix;
    xx=toc;
    fprintf('%03d %f %f %f\n',iter,alpha,NeighSigma,xx);    
end

% All done, and now pack up the results to be sent back to the calling
% function.

results.SOM = SelfOMap;

[results.IDX results.WTS] = SOM_FindClosest(theDataN,SelfOMap);

clear SOM_FindClosest;

results.dataBySOM = theDataN*SelfOMap;

clear theDataN;
clear SelfOMap;

return

%
% All done.
%
