% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
%
%
% Routine to make a map/index of nearest neighbors
%
% function results = SOM_Nearest(nGrid)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function results = SOM_Nearest(nGrid)

% Only look at the relationship's one way - no double counting.

hole = ones(3,3);

hole(2,2) = 0;
hole(2:3,1) = 0;
hole(3,1:2) = 0;

superGrid = zeros(nGrid+2,nGrid+2,nGrid*nGrid);

xx = reshape([1:nGrid^2],[nGrid nGrid]);

for ix = 1:nGrid
    for iy = 1:nGrid
        superGrid(ix:ix+2,iy:iy+2,xx(ix,iy)) = hole;
    end
end

results = squeeze(reshape(superGrid(2:end-1,2:end-1,:),[nGrid nGrid nGrid nGrid]));

%
% alll done
%