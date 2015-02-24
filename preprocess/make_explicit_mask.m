function [] = make_explicit_mask(varargin)

% [] = make_explicit_mask;
%
% If no string input is supplied, runs in the current directory and loads
% the PP.mat file created by preproc_SPM8. Otherwise runs in the directory
% basedir.
%
% Creates a skull-stripped mask from the EPI image and the gray matter
% segmentation
% 
% Optional switches:
% 
% 'erosions' is the number of erosions carried out to attempt to remove
% skull and hair from the EPI. Defaults to 3.
% 
% 'startdir' is the matlab_code directory for the study. Defaults to pwd.
% 
% 'threshP' will modify the automatically determined threshold for removing
% out-of-brain voxels. The threshold will be multiplied by this value, so
% it defaults to one.
% 
% 'wholeimg' will include all clusters, rather than just the largest
% (presumed to be brain). This is useful if the mask is coming out missing
% one or both temporal lobes.
% 
% 'subs' a vector of subject #s to include, in the order they are found in
% the PP.mat for the study.

erosions = 3;
threshP = 1;
wholeimg = 0;

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch(varargin{k})
            case 'erosions'
                erosions = varargin{k+1};
            case 'startdir'
                cd(varargin{k+1})
            case 'threshP'
                threshP = varargin{k+1};
            case 'wholeimg'
                wholeimg = 1;
            case 'subs'
                subs = varargin{k+1};
        end
    end
end

try
    load PP
catch
    try
        cd matlab_code
        load PP
    catch
        error('Can not find PP.mat')
    end
end

if exist('subs')~=1
    subs = 1:size(PP.images,1);
end


write_seg{1}.spm.spatial.normalise.write.roptions.preserve = 0;
write_seg{1}.spm.spatial.normalise.write.roptions.bb = [-78 -120 -80; 78 95 140];
write_seg{1}.spm.spatial.normalise.write.roptions.vox = [3 3 3];
write_seg{1}.spm.spatial.normalise.write.roptions.interp = 1;
write_seg{1}.spm.spatial.normalise.write.roptions.wrap = [0 1 0];
write_seg{1}.spm.spatial.normalise.write.roptions.prefix = 'w3';

smooth_job{1}.spm.spatial.smooth.fwhm = [6 6 6];
smooth_job{1}.spm.spatial.smooth.dtype = 0;
smooth_job{1}.spm.spatial.smooth.im = 0;
smooth_job{1}.spm.spatial.smooth.prefix = 's';

for i = subs
    
    thresh = fmri_mask_thresh_canlab(PP.mean_func{i},'temp.img');
    thresh = thresh * threshP;
    
    [volInfo, dat] = iimg_read_img(PP.mean_func{i},1);
    
    cl = iimg_indx2clusters(dat,volInfo,thresh);
    cl = cl_subdivide(cl,26,1,'corner',1,erosions);
    m = zeros(size(cl));
    if wholeimg
        [imgvec, maskvec] = iimg_clusters2indx(cl,volInfo);
    else
        for j = 1:length(cl)
            m(j) = cl(j).numVox;
        end
        [imgvec, maskvec] = iimg_clusters2indx(cl(m==max(m)),volInfo);
    end
    
    iimg_write_images(imgvec, volInfo, ['..' filesep PP.dirs(i).name filesep 'EPIs' filesep 'functional_mask.img']);
    PP.fmask(i) = filenames(['..' filesep PP.dirs(i).name filesep 'EPIs' filesep 'functional_mask.img'],'absolute');
    
    %     [volInfo, imgvec] = iimg_read_img(fmask{i},1);
    
    write_seg{1}.spm.spatial.normalise.write.subj.matname = {strrep(PP.structurals{i},'.img','_seg_sn.mat')};
    a = strfind(PP.structurals{i},filesep);
    write_seg{1}.spm.spatial.normalise.write.subj.resample = {[PP.structurals{i}(1:a(end)) 'c1' PP.structurals{i}(a(end)+1:end)]};
    spm_jobman('run',{write_seg});
    
    smooth_job{1}.spm.spatial.smooth.data = {[PP.structurals{i}(1:a(end)) 'w3c1' PP.structurals{i}(a(end)+1:end)]};
    spm_jobman('run',{smooth_job});
    
    [v, dat] = iimg_read_img([PP.structurals{i}(1:a(end)) 'sw3c1' PP.structurals{i}(a(end)+1:end)],1);
    iimg_write_images(imgvec > 0.1 & dat > 0.1, volInfo, ['..' filesep PP.dirs(i).name filesep 'stats_mask.img']);
    smask(i) = filenames(['..' filesep PP.dirs(i).name filesep 'stats_mask.img'],'absolute');
end

for i = 1:length(PP.dirs)
    PP.smask{i} = [PP.basedir filesep PP.dirs(i).name filesep 'stats_mask.img'];
    PP.fmask{i} = [PP.basedir filesep PP.dirs(i).name filesep 'EPIs' filesep 'functional_mask.img'];
    if ~exist(PP.smask{i})
        PP.smask{i} = [];
    end
    if ~exist(PP.fmask{i})
        PP.fmask{i} = [];
    end
end
    
    
save PP PP

end

