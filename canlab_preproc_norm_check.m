function CHK = canlab_preproc_norm_check(basedir, varargin)
% CHK = canlab_preproc_norm_check(basedir, subjectwildcard, [anatfilename], [mean_func_filename])
%
% Runs normalization check in standard directories in CANlab preproc
% structure
% Start in directory above individual subjects' data folders.
% T1 images should be in [subjectname]/Structural/SPGR/[anatfilename].img
% mean functional images should be in [subjectname]/Functional/Preprocessed/[meanfuncfilename].img
%
% Subject wildcard should be a cell array of wildcards for directory
% matches. See example:
%
% To use standard names: 'wspgr.img' and 'wmean_ra_func.img'
% CHK = canlab_preproc_norm_check(pwd, {'nsf*' 'NSF*'});
%
% To enter custom names, but use all dirs in current dir for subject list:
% CHK = canlab_preproc_norm_check(pwd, [], 'wspgr.img', 'wmean_ra_func.img');
%
% Tor Wager
% See scnlab_norm_check.m

anatfilename = 'wspgr.img';
mean_func_filename = 'wmean_ra_func.img';
subjectwildcard = '*';

if length(varargin) > 0 && ~isempty(varargin{1})
    subjectwildcard = varargin{1};
end
if length(varargin) > 1 && ~isempty(varargin{2})
    anatfilename = varargin{2};
end

if length(varargin) > 2 && ~isempty(varargin{3})
    mean_func_filename = varargin{3};
end


% Start in directory above individual subjects' data folders.
% T1 images should be in [subjectname]/Structural/SPGR/[anatfilename].img

subjects = canlab_list_subjects(basedir, subjectwildcard{:});

% Standard CANlab file name structure:
[anat, isfile, ~, hasduplicates] = canlab_list_files(subjects, 'Structural', 'SPGR', anatfilename);

if any(~isfile)
    disp('These files do not exist:')
    disp(char(anat{~isfile}))
    disp('Exiting.');
    return
end

if hasduplicates
    disp('Warning!!! Anatomical files have duplicate filenames!');
end

anat = char(anat{:});

mask = which('brainmask.nii');

if ~exist(mask, 'file')
    disp('The mask files does not exist or is not on the path:')
    disp(char(mask{~isfile}))
    disp('Exiting.');
    return
end

template = which('avg152T1.nii');

if ~exist(template, 'file')
    disp('The template file does not exist or is not on the path:')
    disp(char(template{~isfile}))
    disp('Exiting.');
    return
end


% Look for functional mean files:

[func, isfile] = canlab_list_files(basedir, subjects, 'Functional', 'Preprocessed', mean_func_filename);

% If ~exist, try unzipping:
if any(~isfile)
    disp('Cannot find func files. Maybe they are zipped? Trying to unzip.');
    
    for i = 1:length(func)
        eval(['!gunzip ' func{i}]);
    end
    
end

[func, isfile, ~, hasduplicates] = canlab_list_files(basedir, subjects, 'Functional', 'Preprocessed', mean_func_filename);


if any(~isfile)
    disp('These files do not exist:')
    disp(char(func{~isfile}))
    disp('Exiting.');
    return
end

if hasduplicates
    disp('Warning!!! Anatomical files have duplicate filenames!');
end

func = char(func{:});

%%

CHK = scnlab_norm_check(template, anat, func, subjects);

end % function
