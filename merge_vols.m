% merged_vol = merge_vols(files_to_merge, output_base, output_type, [display_commands])
%   Merges a series of Nifti/Analyze volumes into a single 4D volume
%
%   Inputs:
%      files_to_merge - cellstr or char array of filenames
%      output_base - basename of new file (no suffix)
%      output_type - 'NIFTI' or 'ANALYZE' - by default, uses type of first file to be merged
%      display_commands - boolean to display the shell commands (default: 0)
%
%   E.g.:
%      merged_vol = merge_vols(filenames('swravol*nii'), 'swravols', 1);
%      merged_vol = merge_vols(filenames('avol*img'), 'avols', 1);

function merged_vol = merge_vols(files_to_merge, output_base, varargin)
    scn_setup();
    display_commands = 0;

    for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case 'output_type'
                    output_type = varargin{i+1};
                case 'display_commands'
                    display_commands = varargin{i+1};
            end
        end
    end
    if(~exist('output_type', 'var') || isempty(output_type))
        output_type = determine_output_type(files_to_merge{1});
    end
    files_to_merge = cellstr(files_to_merge);
    
    merge_w_fsl(files_to_merge, output_base, output_type, display_commands);
    merged_vol = [output_base determine_ext(output_type)];
end

function output_type = determine_output_type(file_to_merge)
    [pathstr filename ext] = fileparts(file_to_merge);
    switch(ext)
        case {'.img' '.hdr'}
            output_type = 'ANALYZE';
        case '.nii'
            output_type = 'NIFTI';
        otherwise
            error('Unknown extension: %s\n', ext);
    end
end

function ext = determine_ext(output_type)
    switch(output_type)
        case 'ANALYZE'
            ext = '.img';
        case 'NIFTI'
            ext = '.nii';
        otherwise
            error('Unknown output_type: %s\n', output_type);
    end
end

function merge_w_fsl(files_to_merge, output_base, output_type, display_commands)
    global FSLDIR;
    
    merge_prog = sprintf('%s/bin/fslmerge', FSLDIR);
    if(~exist(merge_prog, 'file'))
        fprintf('Could not find fslmerge at %s. Trying avwmerge instead.\n', merge_prog);
        merge_prog = sprintf('%s/bin/avwmerge', FSLDIR);
    end

    merge_command = sprintf('export FSLDIR=%s && . %s/etc/fslconf/fsl.sh && FSLOUTPUTTYPE=%s && %s -t %s %s', ...
        FSLDIR, FSLDIR, output_type, merge_prog, output_base, implode(files_to_merge, ' '));
    if(display_commands), disp(merge_command); end
    unix(merge_command);
end

% function merge_w_afni(files_to_merge, output_base, display_commands)
%     global AFNIDIR;
%
%     [pathstr, name, ext] = fileparts(files_to_merge{1});
%     if(~strcmp(ext, '.hdr'))
%         for i=1:length(files_to_merge)
%             [pathstr, name] = fileparts(files_to_merge{i});
%             files_to_merge{i} = fullfile(pathstr, [name '.hdr']);
%         end
%     end
%
%     if(exist('temp+orig.HEAD', 'file'))
%         unix('rm temp*.HEAD temp*.BRIK');
%     end
%
%     merge_command = sprintf('%s/bin/3dTcat -output temp %s', AFNIDIR, implode(files_to_merge, ' '));
%     if(display_commands), disp(merge_command); end
%     unix(merge_command);
%
%     convert_command = sprintf('%s/bin/3dAFNItoANALYZE -4D %s temp+orig', AFNIDIR, output_base);
%     if(display_commands), disp(convert_command); end
%     unix(convert_command);
%
%     cleanup_command = 'rm temp*.HEAD temp*.BRIK';
%     if(display_commands), disp(cleanup_command); end
%     unix(cleanup_command);
% end

