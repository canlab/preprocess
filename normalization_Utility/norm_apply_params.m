function normdat = norm_apply_params(prm,varargin)
    % normdat = norm_apply_params(prm,['data',dat],['mat',mat],['basis',z,X,Y,BX,BY,BZ,def_flags,d])
    %
    % examples:
    % -----------------------------------------------
    % % This shows you all the steps for applying a warping:
    % name = 'T1.img';
    % V = spm_vol(name);
    % dat = spm_read_vols(V);
    % prm = load('T1_sn.mat');
    % normdat = norm_apply_params(prm,'data',dat)
    %
    % % Write output image:
    % Vout = prm.VG;
    % Vout.fname = 'normalized_t1.img';
    % spm_write_vol(Vout,normdat);
    %
    % Check against template:
    % template = which('avg152T1.mnc');
    % spm_check_registration(str2mat(template,Vout.fname));
    %
    % -----------------------------------------------
    % % This example shows you how you might create the basis set and apply
    % % the warping in separate steps, e.g., for use in another algorithm
    %
    % [z,X,Y,BX,BY,BZ,def_flags,d] = norm_create_basis(prm);
    % normdat = norm_apply_params(prm,'data',dat,'basis',z,X,Y,BX,BY,BZ,def_flags,d);
    %
    % This example shows you how to apply params to an image OTHER THAN the one
    % you used to determine params, by passing in an extra mat matrix:
    % prms = '/Volumes/SCNAlpha/Data_and_Tools/IMAGING_DATA/Opioid_imaging_placebo/Opioid_Placebo3/anatomy/nowarp/spm2-avg152/1001mr_sn.mat';
    % prm = load(prms);
    % name1 = 'mean1001.img';  % your functional to apply to
    % anat = '/Volumes/SCNAlpha/Data_and_Tools/IMAGING_DATA/Opioid_imaging_placebo/Opioid_Placebo3/anatomy/nowarp/spm2-avg152/1001mr.img'
    % V = spm_vol(name1);
    % mymat = V.mat;
    % normdat = norm_apply_params(prm,'data',spm_read_vols(spm_vol(name1)),'mat',mymat);


    % -----------------------------------------------------
    % Process input arguments
    % -----------------------------------------------------
    createBasisFlag = 1;

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'data', dat = varargin{i+1};
                case 'mat', mat = varargin{i+1};
                case 'basis'
                    z = varargin{i+1};
                    X = varargin{i+2};
                    Y = varargin{i+3};
                    BX = varargin{i+4};
                    BY = varargin{i+5};
                    BZ = varargin{i+6};
                    def_flags = varargin{i+7};
                    d = varargin{i+8};
                    createBasisFlag = 0;
                otherwise
                    error('Unknown option');
            end
        end
    end

    % set up image data if not entered
    % -----------------------------------------------------
    if ~exist('dat','var')

        % Read data: map image with nearest neighbor
        % % d  = [def_flags.interp*[1 1 1]' def_flags.wrap(:)];
        % % C = spm_bsplinc(prm.VF,d);
        dat = spm_read_vols(prm.VF);  % seems the same but faster for nearest neighbor
    end

    % Set up affine matrix of image, if not entered
    % -----------------------------------------------------
    if ~exist('mat','var')
        mat = prm.VF.mat;
    end



    % nonlinear basis set: Create (if not input)
    % -----------------------------------------------------
    if createBasisFlag

        [z,X,Y,BX,BY,BZ,def_flags,d] = norm_create_basis(prm);


    end

    Tr = prm.Tr;


    % final affine matrix
    % mat is affine mtx of image
    % VF.mat is affine mtx of image from prm (?difference?)
    % Affine is affine mtx of mapping to template space
    affinemat = mat\prm.VF.mat*prm.Affine;
    %affinemat = mat\prm.VF.mat;  % to leave in native image space

    normdat = zeros(size(X, 1), size(X, 2), length(z));
    for j=1:length(z)   % Cycle over planes
        % Nonlinear deformations
        %----------------------------------------------------------------------------
        tx = get_2Dtrans(Tr(:,:,:,1),BZ,j);
        ty = get_2Dtrans(Tr(:,:,:,2),BZ,j);
        tz = get_2Dtrans(Tr(:,:,:,3),BZ,j);
        X1 = X    + BX*tx*BY';
        Y1 = Y    + BX*ty*BY';
        Z1 = z(j) + BX*tz*BY';

        [X2,Y2,Z2]  = mmult(X1,Y1,Z1,affinemat);
        Dat         = spm_bsplins(dat,X2,Y2,Z2,d);
        %dat(msk{j}) = NaN;

        %normdat(:,:,j) = single(Dat);
        normdat(:,:,j) = Dat;
    end
    normdat = single(normdat);
end


% --------------------------------------------------------------------
%
%
% sub-functions
%
%
% --------------------------------------------------------------------

function T2 = get_2Dtrans(T3,B,j)
    d   = [size(T3) 1 1 1];
    tmp = reshape(T3,d(1)*d(2),d(3));
    T2  = reshape(tmp*B(j,:)',d(1),d(2));
end


function [X2,Y2,Z2] = mmult(X1,Y1,Z1,Mult)
    if length(Z1) == 1,
        X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + (Mult(1,3)*Z1 + Mult(1,4));
        Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + (Mult(2,3)*Z1 + Mult(2,4));
        Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + (Mult(3,3)*Z1 + Mult(3,4));
    else
        X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
        Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
        Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
    end
end
