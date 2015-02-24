% nifti_files = dicom2nifti(dcm_files, num_slices_per_vol, output_base, ['disdaqs', num_disdaqs], ['stack4d' | 'stack3d'], ['verbose'])
%
% Read a specified list of DICOM images (one slice per file), stack them into volumes,
% convert to Analyze format and write an Analyze image/header pair for each volume
%
% Optional inputs:
%    stack3d - split all files into multiple 3D volumes - default
%    stack4d - concatenate all files into one 4D volume
%    disdaqs - removes a certain number of volumes from the beginning of the run - defaults to 0
%    verbose - print out commands
%
% Sample usage:
%   %for functional images
%   dcm_files = filenames('IM*.dcm');
%   num_slices_per_vol = 34;
%   num_disdaqs = 4;
%   output_base = 'vol';
%   nifti_files = dicom2nifti(dcm_files, num_slices_per_vol, output_base, 'disdaqs', num_disdaqs);
%
%   %for structural images
%   dcm_files = filenames('IM*.dcm');
%   num_slices_per_vol = length(dcm_files);
%   output_base = 'spgr';
%   nifti_file = dicom2nifti(dcm_files, num_slices_per_vol, output_base, 'stack3d');
%
% Assumes each DICOM file is a single slice
%
% NB: The number of slices per volume must be passed in as a parameter
% since there is no scanner-independent information to that effect in the
% DICOM header.

function nifti_files = dicom2nifti(dcm_files, num_slices_per_vol, output_base, varargin)
    if(~exist('dcm_files', 'var') || isempty(dcm_files))
        error('No DICOM files to convert');
    end
    if(~exist('num_slices_per_vol', 'var') || isempty(num_slices_per_vol))
        error('num_slices_per_vol not entered');
    end
    if(~ischar(output_base))
        error('Unknown output_base');
    end


    stack4d = 0;
    num_disdaqs = 0;
    verbose = 0;

    if(ischar(dcm_files))
        dcm_files = cellstr(dcm_files);
    end

    for i=1:length(varargin)
        if (ischar(varargin{i}))
            switch(varargin{i})
                case 'stack4d'
                    stack4d = 1;
                case 'stack3d'
                    stack4d = 0;
                case 'disdaqs'
                    num_disdaqs = varargin{i+1};
                case 'verbose'
                    verbose = 1;
            end
        end
    end


    if(rem(length(dcm_files), num_slices_per_vol) ~= 0)
        error('The number of files (%d) is not a multiple of the number of slices per volume (%d)', length(dcm_files), num_slices_per_vol);
    end

    dcm_files = remove_disdaqs(dcm_files, num_disdaqs, num_slices_per_vol, verbose);
    nifti_files = convert_to_nifti(stack4d, dcm_files, num_slices_per_vol, output_base, verbose);
end




function remaining_dcm_files = remove_disdaqs(dcm_files, num_disdaqs, num_slices_per_vol, verbose)
    num_slices_to_move = num_disdaqs * num_slices_per_vol;
    if(num_slices_to_move > 0)
        mkdir('disdaqs');
        mv_command = sprintf('mv %s disdaqs/', implode(dcm_files(1:num_slices_to_move), ' '));
        if(verbose), disp(mv_command); end
        system(mv_command);
    end
    fprintf('Moved %d disdaq slices to disdaqs/\n', num_slices_to_move);
    remaining_dcm_files = dcm_files(num_slices_to_move+1:end);
end


function nifti_files = convert_to_nifti(stack4d, dcm_files, num_slices_per_vol, output_base, verbose)
    global MEDCON;

    scn_setup();

    num_files = length(dcm_files);
    num_volumes = num_files / num_slices_per_vol;

    fprintf('%d total slices\n', length(dcm_files));
    fprintf('%d slices per volume\n', num_slices_per_vol);

    if(stack4d)
        fprintf('%d volumes to be written as one 4-D volume\n', num_volumes);

        temp_base = 'temp';
        generate_3d_vols(dcm_files, num_slices_per_vol, temp_base, 'nifti', verbose); %needed because otherwise there are too many args on the command-line :)
        
        %if(exist(fullfile(pwd, [output_base '.nii']), 'file'))
        %    delete(fullfile(pwd, [output_base '.nii']));
        %end
        merge_command = sprintf('%s -w -stack4d -qc -spm -noprefix -w -c nifti -o %s -f %s*.nii', MEDCON, output_base, temp_base);
        if(verbose), disp(merge_command); end
        system(merge_command);
        nifti_files = {[output_base '.nii']};

        cleanup_command = sprintf('rm %s*.nii', temp_base);
        if(verbose), disp(cleanup_command); end
        system(cleanup_command);
    else
        fprintf('%d volumes to be written as separate 3-D volumes\n', num_volumes);
        nifti_files = generate_3d_vols(dcm_files, num_slices_per_vol, output_base, 'nifti', verbose);
    end

    nifti_files = nifti_files(:);
    fprintf('\nDone\n');
end

function nifti_files = generate_3d_vols(dcm_files, num_slices_per_vol, output_base, output_format, verbose)
    global MEDCON;

    num_files = length(dcm_files);
    num_volumes = num_files / num_slices_per_vol;

    start_indices = 1:num_slices_per_vol:num_files;
    if(length(start_indices) ~= num_volumes), error('The number of volumes does not seem to match in dicom2nifti().'); end

    for i=1:num_volumes
        start_nbr = start_indices(i);
        filelist = implode(dcm_files(start_nbr:start_nbr+num_slices_per_vol-1, :), ' ');

        if(num_volumes > 1)
            current_output_base = sprintf('%s%04d', output_base, i);
        else
            current_output_base = output_base;
        end
        %if(exist(fullfile(pwd, [current_output_base '.nii']), 'file'))
        %    delete(fullfile(pwd, [current_output_base '.nii']));
        %end
        stack_command = sprintf('%s -w -stack3d -qc -spm -noprefix -fh -fv -w -c %s -o %s -f %s', MEDCON, output_format, current_output_base, filelist);
        if(verbose), disp(stack_command); end

        status_string = sprintf('Converting...%d/%d', i, num_volumes);
        fprintf(status_string);
        system(stack_command);
        nifti_files{i} = [current_output_base '.nii'];

        if(i ~= num_volumes)
            erase_string(status_string);
        else
            fprintf('\n');
        end
    end
end