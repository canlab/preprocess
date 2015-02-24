% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
%
%
% function results = SOM_FindClosest(theData,theSOM)
%
% IT IS ASSUMED THAT THE DATA AND SOM ARE UNIT NORMED
%
% theData = theData(nSpace,nTime);
% theSOM  = theSOM(nTime,nSOM);
%
% idx     = array of indices for best matching SOM vector.
% wts     = cos(angle)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [idx, wts] = SOM_FindClosest(theData,theSOM)

dataBySOM = theData*theSOM;

idx  = [];
nvec = [];

[y idx] = sort(dataBySOM,2);

idx = idx(:,end);
wts = y(:,end);

clear theData
clear theSOM

return

%
% All done.
%
