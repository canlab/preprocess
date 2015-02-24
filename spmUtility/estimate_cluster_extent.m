function cl_ext = estimate_cluster_extent(corrected_p, prim_p, con_files, varargin) 

% cl_ext = estimate_cluster_extent(corrected_p, prim_p, con_files, varargin) 
%
% Overview: this function is designed to estimate a cluster extent size for 
% the correction for multiple comparisons based on a Gaussian Random Field 
% Theory using SPM toolboxes. 
%
% Inputs:
% - corrected_p = corrected p value
%   e.g.) cluster-extent corrected p < .05: corrected_p = .05
% - prim_p = primary threshold for height (i.e., cluster-defining threshold)
%   e.g.) prim_p = [0.01 0.005 0.001];
% - con_files = contrast image file names; it could be a cell structure or strings 
%   e.g.) con_file ={'...directory/con0001.img', '...directory/con0002.img',...}
% 
% Output: 
% cl_ext - cl_ext consists of two columns; the first column is the cluster 
%       size that makes a corrected p value under corrected_p (e.g., 0.05). 
%       The second column is corrected p value given the primary height
%       threshold level and the cluster size. The number of rows will be
%       the number of primary p-values
%
% Options:
% - mask = this function will estimate a cluster size for the correction for multiple
%   comparisons within the mask. You can put in a ROI mask or gray matter,
%   whatever. If you don't specify a mask image, brainmask.nii will be
%   used, but the image has to be in your path.
%   e.g.) mask = fullfile(basedir, 'ROI_image.img'); 
%   mask = which('scalped_avg152T1_graymatter_smoothed.img'); % limited to gray matter
%
% ! This will make a temporary directory and delete it in the present directory (pwd).
%
% (usage example)
%   cl_ext = estimate_cluster_extent(0.05, [0.01 0.005 0.001], con_files, 'mask', mask)
%
% by Wani (Choong-Wan) Woo, CANlab, 08/13/12
 
outputdir = pwd;
temp_dir = fullfile(outputdir, 'cl_extent_estimation');
mkdir(temp_dir);
cd(temp_dir);
mask = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'mask', mask = varargin{i+1};
        end
    end
end

if isempty(mask)
    try
        mask = which('brainmask.nii');
    catch
        warning('We tried to use a default mask (brainmask.nii), but we couldn''t use it. Maybe you need to add "brainmask.nii" into your path. We are going to use the first contrast image as a mask.');
        mask = conimgs{1};
    end
end

matlabbatch{1}.spm.stats.factorial_design.dir{1} = temp_dir; % result dir
if ~iscell(con_files)
    for i = 1:size(con_files,1)
        conimgs{i} = deblank(con_files(i,:));
    end
    conimgs = conimgs';
else
    conimgs = con_files;
end

con_num = length(conimgs);

matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = conimgs; % each con image for 2nd level
% matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', [], 'cname',[],'iCFI',[],'iCC',[]);

matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;

matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = mask;

matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

save spm_2nd matlabbatch;

ndf = [1 con_num-1];

for i = 1:length(prim_p)
    u(i) = spm_u(prim_p(i), ndf, 'T');
end

cd(temp_dir);
if exist('SPM.mat', 'file'), delete('SPM.mat'); end
spm_jobman('run', 'spm_2nd.mat');

load SPM.mat;

try 
    cl_ext_spm_spm(SPM); % you need cl_ext_spm_spm.m in your path
catch
    error('You need "cl_ext_spm_spm.m" in your path.');
end

mask_img = filenames('beta*.img','char');
Res_imgs = filenames('ResI_*.img', 'char');
[fwhm, dummy, r] = spm_est_smoothness(Res_imgs, mask_img, [con_num con_num-1]);

V2R = 1/prod(fwhm);

cl_ext = zeros(3,2);

scrsz = get(0, 'ScreenSize');
figure('Position', [1 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/2]);

for i = 1:length(prim_p) % cluster-extent threshold (k) based on RF in SPM

    P = [];
        
    for j = 1:10000
        k = j*V2R; % convert the number os voxels into the number of resels
        
        P(end+1,1) = spm_P_RF(1,k,u(i),ndf,'T',r,1);
        P(end,2) = j;

        if P(end,1) <= corrected_p
            cl_ext(i,1) = j;
            cl_ext(i,2) = P(end,1);
            break
        end
        
    end
    
    hh(i) = subplot(1,length(prim_p),i);
    plot(P(:,2), P(:,1), '-b', 'LineWidth', 1.5);
    xlabel('cluster extent size', 'FontSize', 16);
    ylabel('corrected P value', 'FontSize', 16);
    eval(['title(hh(i), ''Primary P:' num2str(prim_p(i)), ''', ''fontsize'', 16);']);
    set(gca, 'FontSize', 15);
    hold on;
    plot(cl_ext(i,1), cl_ext(i,2), 'r+', 'MarkerSize', 15);
    eval(['text(cl_ext(i,1)/2, cl_ext(i,2), ''+ cl size:' num2str(cl_ext(i,1)) ''', ''FontSize'', 14);']);

end
    
rmdir(temp_dir, 's');
cd(outputdir);

return

