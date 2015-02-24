% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
% Ann Arbor MI.
%
% function results = SOM_Alpha(iter,nIter)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function results = SOM_Alpha(iter,nIter);

global SOM

results = SOM.alpha*exp(-iter/nIter/SOM.learningTimeConstant);

return

%
% All done.
%