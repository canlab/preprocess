function cl = cluster_mask(cl,mask,varargin)
%
% Tor Wager
%
% Masks a clusters structure with an anatomical mask of 1s and 0s
% cl is clusters structure
% mask is string filename of mask (see spm_get, e.g., 'my_mask.img')

tol = 3;    % dominance metric tolerance, in mm
doplot = 0;

V = spm_vol(mask); v = spm_read_vols(V);

[x,y,z] = ind2sub(size(v), find(v));    % find(v==0 | isnan(v)));
XYZmm = voxel2mm([x y z]',V.mat)';
fprintf(1,'Matching clusters ')

if doplot
    spm_image('init',mask)
end



for i = 1:length(cl)
    
    fprintf(1,'%3.0f',i);
    
    [tmp,wh] = dominance_point_match(cl(i).XYZmm',XYZmm,tol);
    
    cl(i).XYZ = cl(i).XYZ(:,wh);
    cl(i).XYZmm = cl(i).XYZmm(:,wh);
    cl(i).Z = cl(i).Z(:,wh);
    
    if doplot & ~isempty(cl(i).XYZ)
        V = cl(i); spm_orthviews('AddColouredBlobs',1,V.XYZ,V.Z,V.M,rand(1,3));
        keyboard
    end
end
        
wh = ones(size(cl));
for i = 1:length(cl)
    if isempty(cl(i).XYZ), wh(i) = 0;,end
end

cl = cl(find(wh));

fprintf(1,'\n');


return
