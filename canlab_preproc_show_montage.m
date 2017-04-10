function canlab_preproc_show_montage(imagename, savefilename)
% canlab_preproc_show_montage(imagename, savefilename)
%
% Display a montage of an image
% Specialized format for canlab_preproc to show image quality.
%
% For more flexible tools for displaying results, etc., see
% fmridisplay.m and cluster_* tools.

% spm_image('init', imagename);
% cluster_orthviews_montage(10, 'axial', [], 'onerow');
% fh = findobj('Tag', 'Graphics'); set(fh, 'Visible', 'off');
% fh = findobj('Tag', 'montage_axial');
% figure(fh);
% scn_export_papersetup(300);
% saveas(fh, savefilename);

close all;
canlab_preproc_montage_first_volumes(imagename)

if nargin > 1
    scn_export_papersetup(800);
    h = findobj('Type','figure');
    if numel(h)>1
        [d, f] = fileparts(savefilename);
        for i = 1:numel(h)
            saveas(h(i), fullfile(d, [f num2str(i) '.png']));
        end
    else
        saveas(h, savefilename);
    end
    
end

end