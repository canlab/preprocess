function canlab_organize_data(basedir,funcdir,structdir,varargin)
%usage: function canlab_organize_data(basedir,funcdir,structdir,{'runappend','dcmconvert',[struct_header]})
%author: Scott Schafer
%date: 6/18/2010
%purpose:  This function organizes fmri data and sets it up to be processed
%          via canlab_preproc.
%
%          Inputs:
%          basedir - This is the primary location of the data.  All file
%                    movements will happen within this directory. This
%                    should be a char array.
%          funcdir -  This should be a vertical cell array of run
%                    directories where the functional images are located.
%          structdir - This should be a char array of the current location
%                    of the structural fmri data
%          varargin:
%               'runappend': Use this if your run directories have very
%                       different names so they can all be called with the
%                       same wildcard ('run*') within canlab_preproc
%               'dcmconvert': If your structural images are in .dcm format,
%                       use this to convert them to a single .nii file.  If
%                       you call this, the next element of varargin must be
%                       the new struct_header that you want to name the
%                       file.
%
dcmconvert = 0;
runappend = 0;
numdirs = numel(funcdir);

for n = 1:numel(varargin)
    switch varargin{n}
        case {'dcmconvert' 'dcm'}
            dcmconvert = 1;
            struct_header = varargin{n+1};
        case 'runappend'
            runappend = 1;
    end
end

if runappend
    rundir = cell(numdirs,1);
    for n = 1:numdirs
        rundir{n} = ['run' funcdir{n}];
    end
else
    rundir = funcdir;
end

%first check and make all required directories
if ~exist(basedir,'dir')
    error('The base directory does not not exist.  Where are you keeping your data?');
end

if ~exist(fullfile(basedir,'Functional'),'dir')
    for n = 1:numdirs
        mkdir(fullfile(basedir,'Functional','Raw',rundir{n}));
        mkdir(fullfile(basedir,'Functional','Preprocessed',rundir{n}));
    end
else
    for n = 1:numdirs
        if ~exist(fullfile(basedir,'Functional','Raw',rundir{n}),'dir')
            mkdir(fullfile(basedir,'Functional','Raw',rundir{n}));
        end
        if ~exist(fullfile(basedir,'Functional','Preprocessed',rundir{n}),'dir')
            mkdir(fullfile(basedir,'Functional','Preprocessed',rundir{n}));
        end
    end
end

if ~exist(fullfile(basedir,'Structural'),'dir')
    mkdir(fullfile(basedir,'Structural','SPGR'));
else
    if ~exist(fullfile(basedir,'Functional','Raw'),'dir')
        mkdir(fullfile(basedir,'Functional','Raw'));
    end
end

%move files
struct_files = filenames(fullfile(basedir,structdir,'*'),'char','absolute');
for j = 1:size(struct_files,1)
    [dummy,f,e] = fileparts(struct_files(j,:));
    if isempty(e)
        error('You need to manually add the file extension to your DICOMs, otherwise movefile doesn''t work.');
%         e = '.dcm';
    end
    movefile(deblank(struct_files(j,:)),fullfile(basedir,'Structural','SPGR',deblank([deblank(f) e])));
end

func_files = cell(numdirs,1);
for n = 1:numdirs
    func_files{n} = filenames(fullfile(basedir,funcdir{n},'*.*'),'char','absolute');
    for j = 1:size(func_files{n},1)
        [dummy,f,e] = fileparts(func_files{n}(j,:));
        movefile(func_files{n}(j,:),fullfile(basedir,'Functional','Raw',rundir{n},[f e]));
    end
end

if dcmconvert
    if ~exist('struct_header','var')
        error('You didn''t include a header for your new Structural file, so all the .dcm files are unchanged, though they have been moved to the new folder.')
    end
    canlab_dcm_converter(struct_header,fullfile(basedir,'Structural','SPGR','*.*'));
    strfile = filenames(fullfile(pwd,[struct_header '.*']));
    [d,f,e] = fileparts(strfile{1});
    movefile(deblank(fullfile(d,[f e])), deblank(fullfile(basedir,'Structural','SPGR',[f e])));
end

