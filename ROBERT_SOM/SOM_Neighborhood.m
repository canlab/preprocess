% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
%
% function results = SOM_Neighborhood(nGrid)
% 
%

function results = SOM_Neighborhood(nGrid)

distMat = zeros(nGrid,nGrid,nGrid*nGrid);

xs = reshape((1:nGrid)'*ones(1,nGrid),[nGrid*nGrid 1]);
ys = reshape(ones(nGrid,1)*(1:nGrid),[nGrid*nGrid 1]);

for iGrid = 1:nGrid*nGrid
  distMat(:,:,iGrid) = reshape(sqrt((xs-xs(iGrid)).^2+(ys- ...
                                                    ys(iGrid)).^2),[10 ...
                      10]);
end

results = reshape(distMat,[nGrid nGrid nGrid nGrid]);;

return

%
% All done
%
