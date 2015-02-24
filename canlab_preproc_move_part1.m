function PREPROC = canlab_preproc_move_part1(PREPROC, varargin)
% canlab_preproc_move_part1(PREPROC)
%
% Special function used in canlab_preproc
% Deletes a* images
% Moves ra*, wra*, swra* images to Preprocessed directory.
% (If these are not present yet, skips them.)
%
% Checks that files exist and saves names in PREPROC structure
% for later use.
%
% Can be run iteratively to check for new out-of-place/missing files.
%
% example:
% PREPROC = canlab_preproc_move_part1(PREPROC);
%
% Optional arguments:
% 'keepavols' : keep slice-timing corrected files (default is to remove)

print_header1('Moving files to Preprocessed and saving file locations in PREPROC');

basedir = PREPROC.basedir;
nruns = length(PREPROC.images_per_session);
rundirs = PREPROC.rundirs;

keepavols = any(strcmp(varargin, 'keepavols'));

for i = 1:nruns
    if ~exist(fullfile(basedir,'Functional','Preprocessed', rundirs(i).name),'dir')
        mkdir(fullfile(basedir,'Functional','Preprocessed', rundirs(i).name))
    end
    
    % delete a vols
    % ----------------------------------------------------------------------
    if ~keepavols
        deletethese = fullfile(PREPROC.basedir, 'Functional','Raw', rundirs(i).name, 'a*.img');
        delete(deletethese);
        
        deletethese = fullfile(PREPROC.basedir, 'Functional','Raw', rundirs(i).name, 'a*.nii');
        delete(deletethese);
        
        deletethese = fullfile(PREPROC.basedir, 'Functional','Raw', rundirs(i).name, 'a*.hdr');
        delete(deletethese);
        
        deletethese = fullfile(PREPROC.basedir, 'Functional','Raw', rundirs(i).name, 'a*.mat');
        delete(deletethese);
    end
    
    % move ra (and wra, swra if present) vols and mat files
    % ----------------------------------------------------------------------
    
    tomove_wildcard{i} = fullfile(basedir, 'Functional', 'Raw', rundirs(i).name, ['*ra*.*']);
    tomove{i} = filenames(tomove_wildcard{i},'char','absolute');
    
    for n = 1:size(tomove{i},1)
        [dummy, f, e] = fileparts(tomove{i}(n,:));
        movefile(strtrim(tomove{i}(n,:)),fullfile(basedir,'Functional','Preprocessed', rundirs(i).name,[f e]));
    end
    %
    %        mats{i} = filenames(fullfile(basedir,'Functional','Raw',rundirs(i).name,'ra*.mat'),'char','absolute');
    %     for n = 1:size(mats{i},1)
    %         [junk, mtf, mte] = fileparts(mats{i}(n,:));
    %         movefile(strtrim(mats{i}(n,:)),fullfile(basedir,'Functional','Preprocessed',rundirs(i).name,[mtf mte]));
    %     end
    
    
    %txt{i} = filenames(fullfile(basedir, 'Functional','Raw',rundirs(i).name,'*.txt'),'char','absolute');
    
    txt{i} = filenames(fullfile(basedir, 'Functional', 'Raw', rundirs(i).name,'rp*.txt'),'char', 'absolute');
    
    for n = 1:size(txt{i}, 1)
        myfile = deblank(txt{i}(n, :));
        [dummy, txf, txe] = fileparts(myfile);
        movefile(myfile,fullfile(basedir,'Functional','Preprocessed',rundirs(i).name,[txf txe]));
    end
    
end

%% List filenames in structure for late use
% Check that files exist, and save PREPROC structure
% -------------------------------------------------------------------------
[tomove_wildcard, PREPROC.a_func_files, PREPROC.ra_func_files, PREPROC.wra_func_files, PREPROC.swra_func_files, PREPROC.mvmt_param_files, allnames] = deal({});

for i = 1:nruns
    
    tomove_wildcard{i} = fullfile(basedir, 'Functional', 'Preprocessed', rundirs(i).name, ['*.*']);
    allnames{i} = filenames(tomove_wildcard{i},'char','absolute');
    
    for n = 1:size(allnames{i},1)
        
        myfile = deblank(allnames{i}(n, :));
        [dummy, f, e] = fileparts(myfile);
        
        if strcmp(e, '.txt') && strcmp(f(1:2), 'rp')
            % realignment params
            PREPROC.mvmt_param_files{i, 1} = myfile;
            
        elseif strcmp(e, '.img') || strcmp(e, '.nii')
            
            myfile = check_valid_imagename(myfile, 0);
            if isempty(myfile)
                disp('canlab_preproc warning!!! File does not exist, and it should:');
                fprintf('%s\n', myfile);
                continue
            end
            
            switch f(1)
                case 'a'  % avol - should not happen unless you keep avols.
                    PREPROC.a_func_files{i, 1} = myfile;
                case 'r'  % ravol
                    PREPROC.ra_func_files{i, 1} = myfile;
                case 'w'  % wravol
                    PREPROC.wra_func_files{i, 1} = myfile;
                case 's'  % swravol
                    PREPROC.swra_func_files{i, 1} = myfile;
            end % switch
            
        end % if
        
    end % for file
    
end % for runs

%% Mask - implicit_mask_file

if isfield(PREPROC, 'implicit_mask_file') && exist(PREPROC.implicit_mask_file, 'file')
    
    myfile = deblank(PREPROC.implicit_mask_file);
    newdir = fullfile(basedir, 'Functional', 'Preprocessed');
    
    newfilename = move_image_file(myfile, newdir);
    PREPROC.implicit_mask_file = newfilename;
    
    
else
    PREPROC.implicit_mask_file = [];
end

end % Function



function print_header1(str, str2)

s = '======================================================================================';
len = length(s);

disp('======================================================================================');
disp('=                                                                                    =');
fprintf('= %s%s=\n', str, repmat(' ', 1, max([1 length(s) - length(str) - 3])));
if nargin > 1
    fprintf('= %s%s=\n', str2, repmat(' ', 1, max([1 length(s) - length(str2) - 3])));
end
disp('=                                                                                    =');
disp('======================================================================================');


end




function newimagefile = move_image_file(myfile, newdir)
% move an .img or .nii file, with associated .hdr and .mat for .img
% files, to a new directory.
% do not try to copy onto itself, and exit gracefully if fails.

[mydir, txf, txe] = fileparts(myfile);

if ~strcmp(mydir, newdir) % cannot move onto itself
    [success, message] = movefile(myfile, newdir);
    
    if ~success
        fprintf('File move failed:\n%s\n', message);
    end
    
    myfile = fullfile(mydir, [txf '.hdr']);
    newfilename = fullfile(newdir, [txf '.hdr']);
    [success, message] = movefile(myfile, newfilename);
    
    if ~success
        fprintf('File move failed:\n%s\n', message);
    end
    
    myfile = fullfile(mydir, [txf '.mat']);
    newfilename = fullfile(newdir, [txf '.mat']);
    [success, message] = movefile(myfile, newfilename);
    
    if ~success
        fprintf('File move failed:\n%s\n', message);
    end
    
end

% return output, whether same as before or not.
newimagefile = fullfile(newdir, [txf txe]);

end
