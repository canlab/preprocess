function [xyzout,xyzmm,cl] = sphere_coords(XYZ,r,varargin)
% [xyzout,xyzmm,clusters] = sphere_coords(XYZmm,r,varargin)
%
%
% XYZ   is n x 3 list of sphere centers (mm coordinates)
% v    is radius of sphere (in mm)
% [mat] optional spm .mat matrix containing voxel to mm space transform
%
% tor wager, aug 12, 2005


% get matrix from input or standard file

if length(varargin) > 0
    m = varargin{1};
else
      % get canonical SPM mm coords
        P = which('scalped_single_subj_T1.img');
        V = spm_vol(P);
        m = V.mat;
end

% change things into voxel coords
V.M = m;  XYZ = mm2voxel(XYZ,V);
rv = mean(abs(r ./ diag(V.M(1:3,1:3))'));   % voxel coords


% work in voxel space

for i = 1:size(XYZ,1)
    xyzbox = build_box(XYZ(i,:),rv);
    
    % this sphere center
    xyzi = repmat(XYZ(i,:),size(xyzbox,1),1);
    
    % distance
    d = sum((xyzbox - xyzi) .^ 2,2) .^ .5;
    wh = d <= rv;
    
    % add center back to voxels, remove negative or zero indices
    wh(find(any(xyzbox < 1),2),:) = 0;
    
    % output coords within radius
    xyzout{i} = xyzbox(wh,:);
    
end

% mm output, if we have a mat input argument
xyzmm = [];


for i = 1:size(XYZ,1)
    xyzmm{i} = voxel2mm(xyzout{i}',m)';
end



% clusters output, if asked for

cl = [];
if nargout > 2
    
    
    for i = 1:size(XYZ,1)
		cl(i).title = ['Sphere ' num2str(rv) ' voxels around ' num2str(XYZ(i,:))];
		cl(i).threshold = 1;
        
        if exist('m') == 1
            cl(i).voxSize = diag(m(1:3,1:3))'; 
            cl(i).M = m;
            cl(i).XYZmm = xyzmm{i}';
        end
        
        a = xyzout{i}';
        cl(i).name = [cl(i).title '_' num2str(i) '_' mat2str(size(a,2)) '_voxels'];
		cl(i).numVox = size(a,2);
		cl(i).Z = ones(1,size(a,2));
		cl(i).XYZ = a;
    end
end

return
        
        %if isempty(xyzmm), error('No coordinates in mask.'), end


        
        
        


function xyz2 = build_box(XYZ,rv)

% build box (list of XYZ voxel coords)

lim = round([XYZ - rv; XYZ + rv]);
diffs = diff(lim);

xtmp = prod([diffs(2)+1 diffs(3)+1]);
ztmp = prod([diffs(1)+1 diffs(2)+1]);

x = repmat((lim(1,1):lim(2,1))',xtmp,1);

y = []; for i=1:diffs(2)+1, 
    ytmp = repmat(lim(1,2)+i-1,diffs(1)+1,1); y = [y;ytmp];,
end
y = repmat(y,diffs(3)+1,1);

ztmp = repmat(lim(1,3),ztmp,1);
z = [];
for i = 1:diffs(3)+1
    z = [z; ztmp+i-1];
end

xyz2 = [x y z];

return
