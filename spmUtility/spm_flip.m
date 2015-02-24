% flips images, writes them and adjusts their mat-files
% FORMAT spm_flip(V, overwrite)
% V - a vector of structures containing image volume information.
% overwrite - [0|1] flag to indicate whether or not to overwrite - 0 by default
% For explanation of the elements of the structure see spm_vol
%___________________________________________________________________________
%
% If spm_flip is called without argument, the interactive function spm_get
% is used for getting images.
%
% The data of the images is flipped in voxel-x-direction.
% MODIFIED 7.10.01 BY TOR WAGER TO FLIP X AND Y DIMS.
% The new images are written to img-files with a leading "f", unless overwriting is on
% The associated new mat-files are calculated to obtain a flipping in
% world-x-direction (the left-right-direction, if the respective prior image
% is aligned properly).
%___________________________________________________________________________

% Author:       Thomas Kamer
% Last change:  2001/04/12
% Version:      1.0
% Dependencies: SPM99
%
% Usage of code from Karl Friston
%
% Distributed under GNU General Public License (GPL) as published by the

% Free Software Foundation (Version 2 or higher)
%
% Laboratory for Psychiatric Brain Research,
% Department for Psychiatry,
% University of Bonn

% Altered: Matthew Davidson
% Added overwriting capability

function spm_flip(V, overwrite)
    if ~exist('overwrite', 'var') || isempty(overwrite)
        overwrite = 0;
    end
    
    if nargin==0
        % get image names
        P     = spm_get(Inf,'*.img',{'select images for flipping'});

        % start progress bar
        spm_progress_bar('Init',length(P),'flipping','');

        % circle through images
        for i = 1:length(P)

            % flip and write
            flip(spm_vol(P{i}), overwrite);

            % show progress
            spm_progress_bar('Set',i);
        end

        % end progress bar
        spm_progress_bar('Clear')
    else
        % circle through images
        for i = 1:length(V)

            % flip and write
            flip(V(i), overwrite);
        end
    end
end


function flip(Vi, overwrite)
    % flips image, calculates mat-data, writes all

    % get image data
    Y             = spm_read_vols(Vi,0);

    % prepare name and header of image
    Vo            = Vi;
    if(~overwrite)
        [p,f,e,v] = fileparts(Vi.fname);
        Vo.fname  = fullfile(p,['f' f e]);
    end
    Vo.descrip    = [Vo.descrip ' - flipped'];

    if(Vi.mat(1,1) < 0)
        Vo.mat(1,1) = -1 .* Vo.mat(1,1);
        Vo.mat(1,4) = -1 .* Vo.mat(1,4);
    else
        % flip image data
        Y = flipdim(Y,1);

        % calculate adjusted mat-file data
        Vo.mat(1:3,4) = Vo.mat(1:3,4)+Vo.mat(1:3,1)*Vo.dim(1);
        Vo.mat(2:3,1) = -Vo.mat(2:3,1);
        Vo.mat(1,2:4) = -Vo.mat(1,2:4);
    end
    
    % write flipped image
    Vo = spm_write_vol(Vo,Y);
end