% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
% Ann Arbor MI.
%
% function results = SOM_NEighborWeights(distances,iter,nIter)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function results = SOM_NeighborWeights(distances,iter,nIter)

global SOM

sigma = SOM.sigma*exp(-iter/nIter/SOM.sigmaTimeConstant);

results = exp(-distances.^2/2/sigma^2);

return

%
% All done.
%
