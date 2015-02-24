function PREPROC = scnlab_preproc_part1_newepi(PREPROC)
% scnlab_preproc_part1_newepi is deprecated. Please use preproc_part1 instead. For more info, type:
% help preproc_part1

help preproc_part1
error('scnlab_preproc_part1_newepi is deprecated. Please use preproc_part1.')
    
% PREPROC = scnlab_preproc_part1(PREPROC)
% PREPROC is an object containing these fields:
%   hires_dir: the directory of the highest-res structural image DICOM files
%   hires_basename: e.g., "t1spgr"
%   inplane_dir: the dir of the in-plane/overlay structural image DICOM files
%   inplane_basename: e.g., "t1overlay"
%   func_dirs: list of dirs containing the functional image DICOM files
%   num_slices_per_vol: the number of slices per volume of the functional iamges
%       - required as there is no scanner-independent way of determining that from the DICOM files for the functional images
%
% E.g.:
% PREPROC.hires_dir = 'structural'; %should not have to change much
% PREPROC.hires_basename = 'T1'; %should not have to change much
% PREPROC.inplane_dir = 'structural'; %should not have to change much
% PREPROC.inplane_basename = 'T1inplane'; %should not have to change much
% PREPROC.dir_names = {'r1' 'r2' 'r3' 'r4' 'r5' 'r6'}; % adjust as needed for the number of runs of your study
% PREPROC.num_slices_per_vol = 24;
% PREPROC.num_vols_per_run = [195 199 199 187 193 195]; %includes disdaqs, which will be discarded
% PREPROC.num_disdaq_vols = 3;
% subjdirs = filenames('rea*');
%
% for i=1:length(subjdirs)
%   cd(subjdirs{i});
%   scnlab_preproc_part1(PREPROC);
%   cd('..');
% end
%
% NOTE: like scnlab_preproc_part1, but modified to handle Columbia fMRI
% Center dicom format for EPI images as of Mar 2007.
% NOTE: ***dicom stacking part not included yet.***

global FSLDIR;
global ANA4DTO3D;

scn_setup();

% we may want to skull strip only the spgr for normalization

% DO THIS ON ALL
% 1 - convert from DICOM
%convert_dirs(PREPROC);

% DO THIS STUFF ON STRUCTURALS
% 1 - brain extract
% 2 - segment (optional)
% 3 - homogeneity correct

% brain extraction

% name1 = 'anatomy/t1overlay.img';
% name2 = extract_brain(name1);
%
%
% name1 = 'anatomy/t1spgr.img';
% name2 = extract_brain(name1);
% spm_check_registration(str2mat(name1,name2))


% DO THIS STUFF ON FUNCTIONALS
% 1 - acquisition timing correction
% 2 - realignment to 1st volume of 1st run
% set images beforehand , e.g.:
%   d = dir('i*.img'); func_files = {};
%   [func_files{1:length(d),1}] = deal(d.name);

% Store number of volumes in each functional run for use later
% Assumes number of slices per volume doesn't change over funcs
num_runs = length(PREPROC.dir_names);



% not needed; this is for old setup
% % % Convert from Stimulate to Analyze format
% % disp('Conversion from Stimulate to Analyze format started.')
% % for i=1:num_runs
% %     run_dir = fileparts(PREPROC.func_files{i});
% %     sdt2an(fullfile(pwd(), run_dir, 'stimulate.sdt'), run_dir);
% % end
% % disp('Conversion from Stimulate to Analyze format finished.')


% Convert DICOM to ANALYZE
% convert_epi_dicom_to_analyze


% Stack strucutral images



% Concatenate 3-D .img volumes from all runs
% so that we can do slice timing, realignment with FSL

% image_list = implode(PREPROC.func_files, ' ');
merge_command = sprintf('. %s/etc/fslconf/fsl.sh && %s/bin/avwmerge -t vols.img vol*img', FSLDIR, FSLDIR);
disp(merge_command);
unix(merge_command);
disp('Done merging.');

% Correct for different times of acquisition in slices
% add a to filename
disp('Acquisition timing correction started.')
slicetimer_command = sprintf('. %s/etc/fslconf/fsl.sh && %s/bin/slicetimer -i vols.img -o avols.img -r 1.5', FSLDIR, FSLDIR);
disp(slicetimer_command);
unix(slicetimer_command);
disp('Acquisition timing correction finished.')

% Motion correction
% WARNING: this may take a long time
% add r to filename
disp('Motion correction started.')
mcflirt_command = sprintf('. %s/etc/fslconf/fsl.sh && %s/bin/mcflirt -in avols.img -out ravols.img -plots -refvol 1 -rmsrel -rmsabs -plots', FSLDIR, FSLDIR);
disp(mcflirt_command);
unix(mcflirt_command);
disp('Motion correction finished.')

% file and slot into its appropriate scan dir
%volume = zeros(header_info(1).dim(1:3));
%new_header_info = spm_vol(['ravols.img']);

mkdir('temp');
split_command = sprintf('!%s ravols.img temp', ANA4DTO3D);
disp(split_command);
eval(split_command);
hdr_files = dir('temp/ravols*.hdr');
img_files = dir('temp/ravols*.img');
num_files = length(img_files);
files_copied = 0;
for i = 1:num_runs
    mkdir(['r' int2str(i) '/'], 'disdaqs');
    for j = 1:PREPROC.num_vols_per_run(i)
        files_copied = files_copied + 1;

        % NB: movefile() is insanely slow, do not use it
        if(j <= PREPROC.num_disdaq_vols)
            unix(sprintf('mv temp/%s r%d/disdaqs/ravol%03d.hdr', hdr_files(files_copied).name, i, j));
            unix(sprintf('mv temp/%s r%d/disdaqs/ravol%03d.img', img_files(files_copied).name, i, j));
        else
            unix(sprintf('mv temp/%s r%d/ravol%03d.hdr', hdr_files(files_copied).name, i, j));
            unix(sprintf('mv temp/%s r%d/ravol%03d.img', img_files(files_copied).name, i, j));
        end
    end
    disp(['Moved ' int2str(PREPROC.num_vols_per_run(i)) ' files for run ' int2str(i)]);
end
rmdir('temp','s');






% INLINE
% --------------------------------------------------------------


    function convert_epi_dicom_to_analyze

        for j=1:length(runs)
            cd(PREPROC.dir_names{j});

            fprintf(1,'In directory: \n %s\n', pwd);

            if exist('dicom', 'dir')
                fprintf(1,'dicom subdirectory found: Assuming already done. \n \n');

            else
                mkdir('dicom')

                dcm_files = filenames('*.dcm');
                fprintf(1,'Converting %3.0f files \n\n', length(dcm_files));

                % make ANALYZE images
                scn_dicom_convert(dcm_files, PREPROC.num_slices_per_vol);

                % move to 'dicom' directory
                % update code to make nicer later
                !mv *.dcm dicom/

            end

            cd('..');
        end

    end

end











function convert_dirs(PREPROC)
%mklinks2 should be used whenever the slices are not taken sequentially
global MKLINKS;

cwd = pwd();

try
    tic;

    %hires
    cd(fullfile(cwd, PREPROC.hires_dir));
    p = filenames('*MRDC*');
    scn_dicom_convert(p, size(p, 1));
    movefile('vol001.hdr', fullfile(cwd, 'anatomy', [PREPROC.hires_basename '.hdr']));
    movefile('vol001.img', fullfile(cwd, 'anatomy', [PREPROC.hires_basename '.img']));

    %inplane
    cd(fullfile(cwd, PREPROC.inplane_dir));
    eval(['!' MKLINKS ' ./'])
    p = filenames('*.dcm');
    scn_dicom_convert(p, size(p, 1));
    movefile('vol001.hdr', fullfile(cwd, 'anatomy', [PREPROC.inplane_basename '.hdr']));
    movefile('vol001.img', fullfile(cwd, 'anatomy', [PREPROC.inplane_basename '.img']));

    %functionals
    for i=1:length(PREPROC.func_dirs)
        cd(fullfile(cwd, PREPROC.func_dirs{i}));
        p = filenames('*MRDC*');
        scn_dicom_convert(p, PREPROC.num_slices_per_vol);
    end

    elapsed_time = toc;
    fprintf('%d seconds to convert all DICOM files.\n', elapsed_time);

    cd(cwd);
catch
    cd(cwd);
end







end



function fout = segment(sout)
% segmentation
[dd,ff,ee] = fileparts(sout);
fout = fullfile(dd,['e' ff ee]);
s = ['!/Applications/FSL/fsl/bin/fast -t1 -c 3 -n -os -od ' sout(1:end-4) ' ' sout(1:end-4) '.hdr'];
disp(s)
eval(s);
end



function nonlinear_smooth(outputfilename)
% 'Nonlinear' noise reduction

[dd,ff,ee] = fileparts(outputfilename);
V = spm_vol(outputfilename);v = spm_read_vols(V);
nn = std(v(:)) ./ 3;
sout = fullfile(dd,['s' ff ee]);
s = ['!/Applications/FSL/fsl/bin/susan_smooth ' outputfilename ' ' num2str(nn) ' ' sout ' 0 3D 1 0'];
disp(s)
eval(s);
end