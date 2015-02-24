% set_hdr_current_coords(imgname, [extra_image_list])
%
% Sets the origin to the current location of the crosshairs and optionally
% copies that origin to the list of files given in varargin (adjusted for
% different dims.

function set_hdr_current_coords(img, varargin)
    % SPM5 reads 4-D images, so let's just expand the names here
    % First image is the one to view/set on, rest are those to apply origin
    % to
    % NOTE: could return V structure here for SPM5, but don't for SPM2
    % compatibility...re-load below.
    disp('Mapping images and counting volumes.');
    
    viewimg = expand_4d_filenames(deblank(char(img)), 1);
    Vviewimg = spm_vol(viewimg);

    imgs = cellstr(expand_4d_filenames(img));
    
    if(~isempty(varargin))
        imgs = vertcat(imgs, cellstr(expand_4d_filenames(strvcat(varargin{:}))));
    end

    if(length(imgs) > 1)
        fprintf('Additional volumes to apply origin to: %3.0f\n', length(imgs) - 1);
    end
    
    % SPM5 will return Vviewimg for each volume, so we need just the first
    % one...


    spm_image('init', Vviewimg);
    spm_orthviews('Reposition', [0 0 0]);

    s = input('Press return when the cursor is on the origin.');

    mm_cursor_pos = spm_orthviews('Pos');    % in mm from current origin

    M = Vviewimg.mat;
    mm_new_origin = M(1:3, 4) - mm_cursor_pos;
    M(1:3, 4) = mm_new_origin;

    mm_dims = (M * [Vviewimg.dim(1:3) 1]') - (M * [0 0 0 1]');
    mm_dims(4) = 1;
    orig_scaled = abs((M * [0 0 0 1]') ./ mm_dims);
    orig_scaled = (orig_scaled(1:3));

    % abs_mm_new_origin = abs(mm_new_origin);
    % vox_origin = M\[0 0 0 1]';
    % vox_origin = vox_origin(1:3);

    % voxel_dims = diag(M(1:3, 1:3));
    % vox_neworigin = (mm_cursor_pos - neg_mm_origin) ./ voxel_dims;  % in voxels
    % mm_neworigin = vox_neworigin .* voxel_dims;   % mm from corner
    % mm_neworigin = mm_cursor_pos - neg_mm_origin;



    % R = M(1:3, 1:3);
    % M(1:3, 4) = -R*(vox_origin + R\mm_cursor_pos);
    spm_get_space(viewimg, M); % despite the name, this invocation actually *sets* the affine matrix

    % REMAINING IMAGES
    % ------------------------------------------------------------------------
    if length(imgs) > 1
        % More volumes in first image and/or additional images

        new_str = sprintf('Applying to: %d/%d', 0, length(imgs));
        fprintf(new_str);
        
        for i = 2:length(imgs)
            old_str = new_str;
            new_str = sprintf('Applying to: %d/%d', i, length(imgs));
            erase_and_display(old_str, new_str);
            
            Vcurrent = spm_vol(imgs{i});

            func_mm_dims = (Vcurrent(1).mat * [Vcurrent(1).dim(1:3) 1]') - (Vcurrent(1).mat * [0 0 0 1]');

            % for inverse dimensions, reverse the scale factor
            % whneg = func_mm_dims < 0;
            whneg = sign(func_mm_dims(1:3)) ~= sign(mm_dims(1:3));
            scale_factor = orig_scaled;
            scale_factor(whneg) = 1 - orig_scaled(whneg);
            mm_new_origin = func_mm_dims(1:3) .* scale_factor;

            %         [DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(p(i, :));
            %         NEWORIGIN = round(mm_neworigin' ./ abs(VOX));
            %         [s] = spm_hwrite(p(i, :), DIM, VOX, SCALE, TYPE, OFFSET, NEWORIGIN, DESCRIP);

            M = spm_get_space(imgs{i});
            %         M(1:3, 4) = (-1 * sign(diag(M(1:3, 1:3)))) .* mm_new_origin;
            M(1:3, 4) = -1 * mm_new_origin;

            spm_get_space(imgs{i}, M);
        end

        fprintf('\n');
    end