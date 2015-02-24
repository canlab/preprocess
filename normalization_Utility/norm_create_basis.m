function [z,X,Y,BX,BY,BZ,def_flags,d] = norm_create_basis(prm,varargin)
    % [z,X,Y,BX,BY,BZ,def_flags,d] = norm_create_basis(prm,[nbf])
    %
    % tor wager, adapted from code by John Ashburner (SPM2, 2004 mod)
    %
    % prm is parameter set from spm_normalize.m
    % uses volume info for 'template' in prm.VG to get:
    % a bounding box (bb)
    % voxel sizes (vox)
    % x, y, z locations of image (??); returns only z
    % grid of X, Y locations in image (??); (X, Y)
    % and nonlinear basis set (BX, BY, BZ) for warping
    %
    % This is a set of the basis functions needed to apply parameters in
    % prm.Tr to an image and get a warped output volume.
    %
    % d is a matrix specifying interpolation and wrapping instructions
    % interpolation = nearest neighbor, wrapping = no
    %
    % nbf is optional; # basis functions in x, y, z

    % default flags
    def_flags = struct('interp',1,'vox',NaN,'bb',NaN,'wrap',[0 0 0],'preserve',0);

    % bounding box and voxel size (get from template)
    [def_flags.bb, def_flags.vox] = bbvox_from_V(prm.VG(1));
    %[def_flags.bb, def_flags.vox] = bbvox_from_V(prm.VF); % to leave in
    %native image space


    % spm options for interpolation and wrapping
    d  = [def_flags.interp*[1 1 1]' def_flags.wrap(:)];



    % x, y, z mapping
    [x,y,z] = get_xyzmat(prm,def_flags.bb,def_flags.vox);  % [x,y,z,mat] = get_xy...not sure what this mat is; not quite same as VF.mat

    [X,Y] = ndgrid(x,y);

    
    % number of x, y, z basis functions
    if nargin > 1 && ~isempty(varargin{1})
        nbf = varargin{1};
    else
        nbf = [size(prm.Tr,1) size(prm.Tr,2) size(prm.Tr,3)];
    end
    
    BX = spm_dctmtx(prm.VG(1).dim(1),nbf(1),x-1);
    BY = spm_dctmtx(prm.VG(1).dim(2),nbf(2),y-1);
    BZ = spm_dctmtx(prm.VG(1).dim(3),nbf(3),z-1);


    return
    
    
    
    
    

    %_______________________________________________________________________
function [bb,vx] = bbvox_from_V(V)
    vx = sqrt(sum(V.mat(1:3,1:3).^2));
    if det(V.mat(1:3,1:3))<0, vx(1) = -vx(1); end;

    o  = V.mat\[0 0 0 1]';
    o  = o(1:3)';
    bb = [-vx.*(o-1) ; vx.*(V.dim(1:3)-o)];
    return;
    %_______________________________________________________________________


    %_______________________________________________________________________
function [x,y,z,mat] = get_xyzmat(prm,bb,vox)
    % The old voxel size and origin notation is used here.
    % This requires that the position and orientation
    % of the template is transverse.  It would not be
    % straitforward to account for templates that are
    % in different orientations because the basis functions
    % would no longer be seperable.  The seperable basis
    % functions mean that computing the deformation field
    % from the parameters is much faster.

    % bb  = sort(bb);
    % vox = abs(vox);

    msk       = find(vox<0);
    bb        = sort(bb);
    bb(:,msk) = flipud(bb(:,msk));

    % Adjust bounding box slightly - so it rounds to closest voxel.
    bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
    bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
    bb(:,3) = round(bb(:,3)/vox(3))*vox(3);

    M   = prm.VG(1).mat;
    vxg = sqrt(sum(M(1:3,1:3).^2));
    if det(M(1:3,1:3))<0, vxg(1) = -vxg(1); end;
    ogn = M\[0 0 0 1]';
    ogn = ogn(1:3)';

    % Convert range into range of voxels within template image
    x   = (bb(1,1):vox(1):bb(2,1))/vxg(1) + ogn(1);
    y   = (bb(1,2):vox(2):bb(2,2))/vxg(2) + ogn(2);
    z   = (bb(1,3):vox(3):bb(2,3))/vxg(3) + ogn(3);

    og  = -vxg.*ogn;
    of  = -vox.*(round(-bb(1,:)./vox)+1);
    M1  = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
    M2  = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
    mat = prm.VG(1).mat*inv(M1)*M2;

    if (spm_flip_analyze_images & det(mat(1:3,1:3))>0) | (~spm_flip_analyze_images & det(mat(1:3,1:3))<0),
        Flp = [-1 0 0 (length(x)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
        mat = mat*Flp;
        x   = flipud(x(:))';
    end
    return

