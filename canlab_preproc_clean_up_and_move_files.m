function PREPROC = canlab_preproc_clean_up_and_move_files(PREPROC, varargin)
% canlab_preproc_clean_up_and_move_files(PREPROC)
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

%Programmers' Notes:
% 3/20/13 - Luke Chang - Added a test to ignore move command if no file
% exists


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
    % continue gracefully if there aren't any ra files to move and checks Preproc folder - 3/20/13 - Luke
    % ----------------------------------------------------------------------
    raexist = dir(fullfile(rundirs(i).folder, rundirs(i).name,'ra*'));
    
    if ~isempty(raexist)
        
        tomove_wildcard{i} = fullfile(rundirs(i).folder, rundirs(i).name, '*ra*.*');
        tomove{i} = filenames(tomove_wildcard{i},'char','absolute');
        
        for n = 1:size(tomove{i},1)
            [dummy, f, e] = fileparts(deblank(tomove{i}(n, :)));
            movefile(strtrim(tomove{i}(n,:)),fullfile(basedir,'Functional','Preprocessed', rundirs(i).name,[f e]));
        end
        
    else  %check if it already exists in Preprocessed folder.
        
        raexist = dir(fullfile(basedir,'Functional','Preprocessed', rundirs(i).name, 'ra*'));
        
        if isempty(raexist)
            error('RA files do not exist in Raw or Preprocessed folders.  Run Preproc Part 1 again')
        end
        
    end

    % Move alternative task-independent realigned files (TIR_*) if they
    % exist. These are not part of the standard preproc path as of Jan
    % 2014, so they may not exist.
    % move ra (and wra, swra if present) vols and mat files
    % continue gracefully if there aren't any ra files to move and checks Preproc folder - 3/20/13 - Luke
    % ----------------------------------------------------------------------
    TIRexist = dir(fullfile(basedir,'Functional','Raw', rundirs(i).name,'TIR*'));
    
    if ~isempty(TIRexist)
        
        tomove_wildcard{i} = fullfile(rundirs(i).folder, rundirs(i).name, ['*TIR_*.*']);
        tomove{i} = filenames(tomove_wildcard{i},'char','absolute');
        
        for n = 1:size(tomove{i},1)
            [dummy, f, e] = fileparts(deblank(tomove{i}(n, :)));
            movefile(strtrim(tomove{i}(n,:)),fullfile(basedir,'Functional','Preprocessed', rundirs(i).name,[f e]));
        end
        
    end
    
    
    % Realignment params
    % ----------------------------------------------------------------------

    %txt{i} = filenames(fullfile(basedir, 'Functional','Raw',rundirs(i).name,'*.txt'),'char','absolute');
    
    rpexist = dir(fullfile(rundirs(i).folder, rundirs(i).name,'rp*.txt'));
    
    if ~isempty(rpexist)
        
        txt{i} = filenames(fullfile(rundirs(i).folder, rundirs(i).name,'rp*.txt'),'char', 'absolute');
        
        for n = 1:size(txt{i}, 1)
            myfile = deblank(txt{i}(n, :));
            [dummy, txf, txe] = fileparts(myfile);
            movefile(myfile,fullfile(basedir,'Functional','Preprocessed',rundirs(i).name,[txf txe]));
        end
        
    else  %check if it already exists in Preprocessed folder.
        
        rpexist = dir(fullfile(basedir,'Functional','Preprocessed', rundirs(i).name,'rp*txt'));
        
        if isempty(rpexist)
            error('RP text file does not exist in Raw or Preprocessed folders.  Run Preproc Part 1 again')
        end
    end
    
     % Move alternative task-independent realigned parame files (rpTIR_*) if they
    % exist. These are not part of the standard preproc path as of Jan
    % 2014, so they may not exist.
    % ----------------------------------------------------------------------
    rpTIRexist = dir(fullfile(basedir,'Functional','Raw', rundirs(i).name,'rpTIR*.txt'));
    
    if ~isempty(rpTIRexist)
        
        txt{i} = filenames(fullfile(rundirs(i).folder, rundirs(i).name,'rpTIR*.txt'),'char', 'absolute');
        
        for n = 1:size(txt{i}, 1)
            myfile = deblank(txt{i}(n, :));
            [dummy, txf, txe] = fileparts(myfile);
            movefile(myfile,fullfile(basedir,'Functional','Preprocessed',rundirs(i).name,[txf txe]));
        end
        
    end
    
end  % runs


%% List filenames in structure for later use
% Check that files exist, and save PREPROC structure
% -------------------------------------------------------------------------
[tomove_wildcard, PREPROC.a_func_files, PREPROC.ra_func_files, PREPROC.wra_func_files, PREPROC.swra_func_files, PREPROC.mvmt_param_files, allnames] = deal({});

for i = 1:nruns
    
    tomove_wildcard{i} = fullfile(basedir, 'Functional', 'Preprocessed', rundirs(i).name, ['*.*']);
    allnames{i} = filenames(tomove_wildcard{i},'char','absolute');
    
    for n = 1:size(allnames{i},1)
        
        myfile = deblank(allnames{i}(n, :));
        [dummy, f, e] = fileparts(myfile);
        
        if strcmp(e, '.txt') && strcmp(f(1:2), 'rp') && ~strcmp(f(1:5), 'rpTIR')
            % realignment params
            PREPROC.mvmt_param_files{i, 1} = myfile;

        elseif strcmp(e, '.txt') && strcmp(f(1:5), 'rpTIR')
            % realignment params - TIR
            PREPROC.TIR_mvmt_param_files{i, 1} = myfile;
            
        elseif strcmp(e, '.img') || strcmp(e, '.nii')
            
            myfile2 = check_valid_imagename(myfile, 0);
            if isempty(myfile2)
                disp('canlab_preproc warning!!! File does not exist, and it should:');
                fprintf('%s\n', myfile);
                continue
            else
                myfile = myfile2;
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
                case 'T'  % TIR - task-independent realignment file
                    PREPROC.TIR_func_files{i, 1} = myfile;
                    
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
    
    myfile = fullfile(mydir, txf, '.hdr');
    newfilename = fullfile(newdir, [txf '.hdr']);
    [success, message] = movefile(myfile, newfilename);
    
    if ~success
        fprintf('File move failed:\n%s\n', message);
    end
    
    myfile = fullfile(mydir, txf, '.mat');
    newfilename = fullfile(newdir, [txf '.mat']);
    [success, message] = movefile(myfile, newfilename);
    
    if ~success
        fprintf('File move failed:\n%s\n', message);
    end
    
end

% return output, whether same as before or not.
newimagefile = fullfile(newdir, [txf txe]);

end
