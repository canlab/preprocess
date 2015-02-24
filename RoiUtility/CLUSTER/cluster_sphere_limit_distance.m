function XYZ = cluster_sphere_limit_distance(XYZ,d,varargin)
% function XYZ = cluster_sphere_limit_distance(XYZ,d,[opt] MAT,[opt] Z)
%
% given XYZ coordinate list (in 3 row vectors),
% this function returns a list of local maxima that
% are at least d units apart, where units are units
% of input coordinates.
%
% searches sequentially through the space, eliminating
% local maxima within the specified distance
%
% if a 3rd argument is entered, it should be the mat file info
% from SPM - e.g., what's in VOL.M.
% this indicates that XYZ should be translated to mm before selection.
%
% the 4th argument is Z scores for the XYZ values
% indicating that preference is given to the subclusters with 
% the highest Z scores when 2 or more clusters are within d
% WITHOUT this argument, the 1st coord is always saved, introducing
% possible BIAS.
%
% Tor Wager, 1/30/02

XYZin = XYZ;

% if 3rd argument, then we need to transform to mm, so do it.
if nargin > 2, 
    M = varargin{1};,
    XYZ = voxel2mm(XYZ,M);
end

if nargin > 3,
    Z = varargin{2};
end

Q       = ones(1,size(XYZ,2));
omit    = 0*Q;
    
for i = 1:size(XYZ,2)
    
    if ~omit(i) % skip to the next non-omitted coordinate
        
        j       = sum((XYZ - XYZ(:,i)*Q).^2) <= d^2;
    
        %j(1:i)  = 0; old way
    
        % do not omit the 1st or the max Z coordinate (mark 0)
        if exist('Z') == 1
            mymax = find(Z == max(Z(j))); mymax = mymax(1);
            j(mymax) = 0;
        else
            j(i) = 0;
        end
    
        omit    = omit + j;
    end
    
end

XYZ = XYZin;
XYZ(:,find(omit)) = [];

return