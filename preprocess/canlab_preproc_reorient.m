function canlab_preproc_reorient(imgs)
% canlab_preproc_reorient(imgs)
%
% This is used in canlab_preproc_2012
% It shows you the first image in a functional series,
% lets you reorient it in SPM, and then applies the transformation to all
% the images in the series.
%
% Imgs: A cell array, one cell per run, with string matrices in each one.
% Images must be expanded filenames (see expand_4d_filenames.m) with 3-D
% volume names.
%
% Used to check orientation of functional files in PREPROC.func_files
% before other processing stages.
%
% Tor Wager, Aug 2012

%mat = spm_matrix(st.B);

global st

% Show first volumes:
canlab_preproc_montage_first_volumes(imgs)

% Bring up first image:
spm_image('init', imgs{1}(1, :));

drawnow
try_snapnow_for_publish

disp('Check images are in correct orientation.')
disp('Coronal should be top left, saggital should be top right, axial bottom left.')
disp('If there is a marker, it should be on the correct side (usually R in canlab.)');
disp('The saggital slice should be oriented with the eyes to the left,');
disp('and the axial with the eyes towards the top');

% SPECIAL non-interactive flip of all:
% st.B = [0 0 0 0 0 0 -1 1 1 0 0 0];

input('Reorient image and press return to apply to *all* images entered: ');


if ~any(st.B - [0 0 0 0 0 0 1 1 1 0 0 0])
    disp('No reorientation requested. Proceeding without altering the images.');
    return
end

%% Read existing orientations

% in SPM orthviews, st.B contains the orientation info entered into the GUI
% for reorientation purposes.

mat = spm_matrix(st.B);
%mat = spm_matrix([0 0 0 0 0 0 1 1 1 0 0 0]);

P = char(imgs{:});

Mats = zeros(4,4,size(P,1));

for i=1:size(P,1)
    Mats(:,:,i) = spm_get_space(P(i,:));
end

%% Set new orientation
%%

for i=1:size(P,1)
    spm_get_space(P(i,:),mat*Mats(:,:,i));
end

% Show us again.
canlab_preproc_montage_first_volumes(imgs)

% Bring up first image:
spm_image('init', imgs{1}(1, :));
drawnow
try_snapnow_for_publish
    
end % main function




function try_snapnow_for_publish

fhi = findobj('Type', 'Figure', 'Tag', 'Interactive');
if ishandle(fhi), set(fhi, 'Visible', 'off'); end

try
    snapnow
catch
    warning('Error executing snapnow.  Old version of Matlab??');
end


end


