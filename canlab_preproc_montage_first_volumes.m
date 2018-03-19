function canlab_preproc_montage_first_volumes(imgs)
% canlab_preproc_montage_first_volumes(imgs)
%
% Show first functional volume
%
% e.g.,
% canlab_preproc_montage_first_volumes(PREPROC.func_files) % before realignment, etc.
%
% 09/01/12 modified by Wani to make images bigger when the session number 
%          is bigger than 5. 

dosplit = false;

if ~iscell(imgs)
    imgs = cellstr(imgs);
end

n_imgs = length(imgs);

if n_imgs > 6
    dosplit = true;
    n = 6;
else
    n = n_imgs;
end

for i = 1:n_imgs
    imgs{i} = deblank(imgs{i}(1, :));
    
    wh = find(imgs{i} == ',');
    
    if isempty(wh)
    imgs{i} = [deblank(imgs{i}) ',1']; % first volume
    end
end

create_figure('First_volumes', n, 1);

% enforce one row, and set scaling across all images
% use mask that is one image, so hopefully show all voxels
% Wani: sometimes, imgs are in a different space. 
try 
    dat = fmri_data(char(imgs{:}), imgs{1}); 
    d = dat.dat(:); 
    spacing = ceil(dat.volInfo.dim(3) ./ 12);
    spacing = repmat(spacing, n, 1);
catch
    dat_dat = [];
    for i = 1:n
        temp = fmri_data(char(imgs{i}), imgs{i}); 
        dat_dat = [dat_dat; temp.dat];
        spacing(i) = ceil(temp.volInfo.dim(3) ./ 12);
    end
    d = dat_dat(:); 
end

d(d == 0 | isnan(d)) = [];
clim = [prctile(d, 5) prctile(d, 95)];

for i = 1:n
    
    dat = fmri_data(imgs{i}, imgs{i});
    
    axh = subplot(n, 1, i);
    
    display_slices(dat, 'axial', 'spacing', spacing(i), 'slices_per_row', 12, 'voxelspace'); % every 4th slice
    
    title(['Run ' num2str(i)])
    
    set(axh, 'CLim', clim);
    
    drawnow
end

sz = get(0, 'screensize');
set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85*n/6]); % Wani added n/6 to make this proportional to n.

if dosplit
    create_figure('First_volumes_cont', n_imgs-n, 1);
    
    % enforce one row, and set scaling across all images
    % use mask that is one image, so hopefully show all voxels
    
    for i = n+1:n_imgs
        
        dat = fmri_data(imgs{i}, imgs{i});
        
        axh = subplot(n_imgs-n, 1, i-n);
        
        display_slices(dat, 'axial', 'spacing', spacing, 'slices_per_row', 12, 'voxelspace'); % every 4th slice
        
        title(['Run ' num2str(i)])
        
        set(axh, 'CLim', clim);
        
        drawnow
    end
    
    sz = get(0, 'screensize');
    set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85*(n_imgs-n)/6]);
    
end
end % function

