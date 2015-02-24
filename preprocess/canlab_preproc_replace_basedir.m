function canlab_preproc_replace_basedir(new_basedir, varargin)
% Replace base directory path names in all PREPROC_SETUP.mat files to set paths on new local computer.
%
% This should be a ONE-TIME RUN script on your local lab computer
% It replaces the path names saved in PREPROC_SETUP.mat files with ones for
% your local computer.
%
% This gives you convenient access to the names of all images saved in
% PREPROC_SETUP files, and allows you to run model scripts on any machine.
%
% new_basedir is the directory ABOVE individual subject directories where
% the canlab_preproc pipeline has been run.  Files must be in standard
% canlab_preproc organization.
% This function operates on all
% [subjectdir]/Functional/Preprocessed/PREPROC_SETUP.mat files and
% REPLACES AND SAVES THE FILES.
%
% Optional inputs:
% canlab_preproc_replace_basedir(new_basedir, ['test']) 
%    -> Do not save files, test mode only.
%
%
% Examples:
% ------------------------------------------------------------------------------
% new_basedir = '/Users/Shared/psyc7215/Data/Data_Wager_NSF_Pain/SubjectData';
% canlab_preproc_replace_basedir(new_basedir)
%
% OR
% Base data dir on Tor's laptop (customize for your own laptop)
% new_basedir = '/Users/tor/Dropbox/psyc7215_class_files/Part2_Machine_Learning/Data/Data_Wager_NSF_Pain/SubjectData';
% canlab_preproc_replace_basedir(new_basedir)

oldpwd = pwd;
cd(new_basedir);

dosave = 1;
if any(strcmp(varargin, 'test')), dosave = 0; end

% Get subject list
subjects = canlab_list_subjects(new_basedir, '*');

% Get file list
[files, exists] = canlab_list_files(subjects, 'Functional', 'Preprocessed', 'PREPROC_SETUP.mat');
files = files(exists > 0);

subjects = subjects(exists > 0);

if isempty(files), disp('NO MATCHING FILES'); return, end

z = '-----------------------------------';

%% Run for each file

for s = 1:length(files)
    
    fprintf('%s\nSubject %s\n%s\n', z, subjects{s}, z)
    
    
    clear PREPROC
    load(files{s})
    
    old_basedir = fileparts(PREPROC.basedir);  % base only, not subject
    
    N = fieldnames(PREPROC);
    
    % Exclude some fields
    omit_fields = {'TR' 'Timeseries' 'Nuisance' 'num_disdaqs' 'run_wildcard' 'image_wildcard'};
    for k = 1:length(omit_fields)
        wh = strfind(N, omit_fields{k});
        for j = 1:length(wh), if isempty(wh{j}), whomit(j, 1) = false; else, whomit(j, 1) = true; end, end
        N(whomit) = [];
    end
    
    for k = 1:length(N)
        myfield = PREPROC.(N{k});
        if isempty(myfield) || (~ischar(myfield) && ~iscell(myfield)) || (iscell(myfield) && ~ischar(myfield{1}))
            whomit(k, 1) = true;
        else
            whomit(k, 1) = false;
        end
    end
    N(whomit) = [];
    
    % Replace strings
    for i = 1:length(N)
        PREPROC.(N{i}) = strrep(PREPROC.(N{i}), old_basedir, new_basedir);
    end
    
    % func_files
    for i = 1:length(PREPROC.func_files)
        for j = 1:length(PREPROC.func_files{i})
            tmp{j, 1} = deblank(strrep(PREPROC.func_files{i}(j, :), old_basedir, new_basedir));
        end
        PREPROC.func_files{i} = char(tmp{:});
        
        fprintf('Run %3.0f\t', i)
        check_valid_imagename(PREPROC.func_files{1}, 2);
        
    end
    
    fprintf('\n');
    
    if dosave 
        save(files{s}, '-append', 'PREPROC');
    end
    
end

cd(oldpwd)

end % function


