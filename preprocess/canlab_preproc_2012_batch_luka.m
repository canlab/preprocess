function canlab_preproc_2012_batch(basedir,runwc,imgwc,anat,TR,acqorder,ndisdaqs,varargin)
%usage: canlab_preproc_2012_batch(basedir,runwc,imgwc,anat,TR,acqorder,ndisdaqs,varargin)
%author: Scott Schafer (adapted from Tor's canlab_preproc_2012_batch_script)
%date: 9/20/2012
%purpose: This is a function that allows a user to preprocess fmri data in 
%         a single batch process in the most efficient way possible.  It
%         assumes that you have the standard canlab directory structure,
%         though it allows for some slight deviations as discussed below.
%
%         There are four main pieces to this script, two interactive and
%         two non-interactive parts.  Both non-interactive parts are
%         published to an html script.
%         Each piece runs for all selected subjects before moving to the next:
%
%           Reorietation (interactive)
%           Quality control, slice timing, motion correction (non-interactive)
%           Origin setting and coregistration (interactive)
%           Warping and smoothing (non-interactive)
%
%         It will also save the following for each subject:
%           - scripts with the overall non-interactive preprocess commands
%             for each subject, to be passed to publish.m, in [datadir]/scripts
%           - A file called PREPROC_SETUP.mat in [subjectdir]/Functional/Preprocessed
%             that contains all the image names, etc. This file can be used directly,
%             and it is loaded in canlab_preproc_2012 before each run, and updated
%             throughout preprocessing. 
%
%         Input variables:
%                basedir - char array of the directory for the main study, 
%                           usually two levels above the individual subjects' 
%                           data directories.  **on dream, may need to pass
%                           in basedir as /dreamio3/wagerlab rather than
%                           /data/projects/wagerlab
%                  runwc - char array that defines a
%                           wildcard to choose runs within the subject
%                           directories, presumably within
%                           [subjdir]/Functional/Raw/
%                           if you have two levels of wildcards in each subj dir, must enter
%                           a  2-row char matrix, i.e. ['time*'; 'run*']
%                  imgwc - char array that defines the wildcard that
%                           selects images within the Raw directory
%                   anat - char array that defines the anatomical image
%                           found in Structural/SPGR/
%                     TR - double value that defines the TR the images were 
%                           acquired at.
%               acqorder - char array defining slice acquisition order.
%                           several options are available, the most common
%                           being 'interleaved_BU' - see canlab_preproc_2012
%                           for more options
%               ndisdaqs - double value that represents the number of
%                           images to throw away prior to analysis.  If
%                           you've already done this, set it to 0!
%         Optional input:
%              'datadir' - follow this flag with a char array of the
%                          directory one level above the individual subject
%                          directories.  By default this is assumed to be
%                          [basedir]/Imaging
%               'subjwc' - follow this flag with a char array that defines
%                          wildcard that identifies subject directories
%                          within datadir.  This is required if
%                          EXPT_batch.mat does not exist within [datadir]. 
%                          EXPT_batch.mat contains an EXPT struct, where
%                          EXPT.subjects is a cell array of subject
%                          directories.
%             'wh_subjs' - follow this flag with a integer array that
%                          defines which subjects should be preprocessed.
%                          If not input, the batch will ask you to define
%                          it based on a list of subjects.
%       'no_interactive' - A flag for not running interactive parts of
%                          canlab_preproc_2012.  This avoids reorientation,
%                          origin setting AND coregistration.
%             'no_part1' - A flag for not running the qc metrics or the
%                          slice timing and motion correction.
%             'no_part2' - A flag for not running the warping or smoothing.
%
%     'no_reorientation' - A flag that allows you to skip reorientation
%                          while still doing origin setting.  You should
%                          always look at your raw images, so only do this
%                          if you've already gone through the reorientation
%                          check at least once.
%
%       'dream'          - If on dream, use this option to run subjects in 
%                          parallel. It's recommended that you do reorientation
%                          between parts 1 and 2. see wiki for more details.
%       'dream_a'        - Same as 'dream' except "automatic" in that it
%                          won't ask for confirmation before submitting
%                          jobs.
%
%       Example usage:
%           canlab_preproc_2012_batch(pwd,'r*','LEVO*','mpr*.nii',1.3,'interleaved_BU',0)
%

%initialize defaults
do_part1 = 1;
do_interactive = 1;
do_part2 = 1;
do_reorient = 1;
dream = 0;
dreamask = 1;

VERBOSE = 0;

%parse varargin
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'datadir'
                datadir = varargin{i+1};
            case 'subjwc'
                subjwc = varargin{i+1};
            case 'wh_subjs'
                subjs_to_run = varargin{i+1};
            case 'no_part1'
                do_part1 = 0;
            case 'no_interactive'
                do_interactive = 0;
            case 'no_part2'
                do_part2 = 0;
            case 'no_reorientation'
                do_reorient = 0;
            case 'dream'
                dream = 1;                
            case 'dream_a'
                dream = 1;
                dreamask = 0;
            otherwise
                warning(sprintf('Unknown character varargin entered: %s  Unknown use, skipping',varargin{i}))
        end
    end
end

%dream flag trumps everything else
if dream
    do_interactive=0;
    do_reorient=0;
end

%check whether datadir entered as a variable
if ~exist('datadir','var')
    datadir = fullfile(basedir,'Imaging');
end

%check whether datadir exists as a directory.  If so, move there
if ~exist(datadir,'dir')
    error('Data directory (%s) does not exist.  Re-check your directory structure.',datadir)
else
    cd(datadir)
end

%check for EXPT_batch that contains information about the subjects in the study
%If it doesn't exist, create it
if exist('EXPT_batch','file')
    fprintf('EXPT_batch file found, loading subjects...\n')
    load EXPT_batch
else
    fprintf('No EXPT_batch file, constructing...\n')
    %Check to see if the user input a wild card to identify subject directories
    if ~exist('subjwc','var')
        error('No subject wildcard! Cannot construct EXPT_batch file')
    else
        %If the subject wildcard is found, construct subject directory list
        subjects = dir(subjwc);
        for i = 1:numel(subjects)
            if subjects(i).isdir
                EXPT.subjects{i,1} = subjects(i).name;
            end
        end
        save EXPT_batch EXPT
    end
end

%check whether the input subject range exceeds dimensions
if exist('subjs_to_run','var')
    if max(subjs_to_run) > numel(EXPT.subjects)
        error('Invalid subject vector - values too large')
    end
    %display which subjects will be preprocessed
    fprintf('The following subjects will be preprocessed:\n')
    for i = 1:numel(subjs_to_run)
        fprintf('%3.0f\t%s\n',i,EXPT.subjects{i})
    end
else
    fprintf('\n')
    %If not specified beforehand, force user to choose which subjects to use.
    for i = 1:numel(EXPT.subjects)
        fprintf('%3.0f\t%s\n',i,EXPT.subjects{i})
    end
    subjs_to_run = input('Please input a vector of subjects you would like to preprocess.\nFor example, [1 2] will run the first two subjects: ');
    
    %Check for invalid ranges
    if max(subjs_to_run) > numel(EXPT.subjects)
        error('Invalid subject vector - values too large')
    end
    
    %display which subjects will be run
    fprintf('The following subjects will be preprocessed:\n')
    for i = subjs_to_run
        fprintf('%3.0f\t%s\n',i,EXPT.subjects{i})
    end
end

%% Begin preprocessing

fprintf('----------------------------------------------------\n')
fprintf('|                                                  |\n')
fprintf('|           Beginning preprocessing                |\n')
fprintf('|                                                  |\n')
fprintf('----------------------------------------------------\n\n\n')

%first, check orientations if doing interactive
if do_interactive && do_reorient
    fprintf('Beginning interactive part1:\n')
    fprintf('Please check all of the following images to make sure they are in the correct orientation\n')
    fprintf('Do not bother setting origins yet, nothing will happen anyway.\n')
    for subj = subjs_to_run
        func_files = canlab_preproc_list_files(fullfile(datadir, EXPT.subjects{subj}), runwc, imgwc, ndisdaqs);
        canlab_preproc_reorient(func_files);
    end
    fprintf('\n')
else
    fprintf('No reorientation check requested\n\n')
end

%setup for publishing
scriptbasedir = fullfile(datadir, 'scripts');
if ~exist(scriptbasedir, 'dir')
    mkdir(scriptbasedir)
end
scriptdir = fullfile(scriptbasedir,datestr(now,'yyyy-mm-dd_HH_MM_SS'));
if exist(scriptdir, 'dir')
    fprintf('deleting existing script directory: %s\n',scriptdir);
    rmdir(scriptdir,'s');
end
mkdir(scriptdir)
addpath(scriptdir)
% if dream, savepath(fullfile(scriptdir,'pathdef.m')); end

%Do the quality control images, slice acquisition timing, and motion correction
if do_part1
    runcommandstring = 'canlab_preproc_2012(subjectdir, runwc, imgwc, TR, anat, ''acquisition'', acqorder, ''nointeractive'',''nopart2'',''nomovie''); \n';
    jobstr = 'preproc_part1';
    fprintf('Initiating preproc_part1: Slice acquisition timing and motion correction\n\n')
    for subj = subjs_to_run
        runcommand(subj,EXPT,scriptdir,jobstr,datadir,basedir,runcommandstring,TR,anat,acqorder,runwc,imgwc,dream);
    end
    fprintf('\nPreproc part 1 completed.\n\n')
else
    fprintf('Skipping preproc part 1\n\n')
end

%Set the origins!
if do_interactive
    fprintf('Set the origins and coregister now\n\n')
    for subj = subjs_to_run
        canlab_preproc_2012(fullfile(datadir,EXPT.subjects{subj}), runwc, imgwc, TR, anat, 'acquisition', acqorder,'nopart1','nowarp','nomean','nosmooth','nospikeid','nomovie','nodatachecks');
    end
else
    fprintf('Skipping origin setting and coregistration...\n\n')
end

%Do smoothing and warping
if do_part2
    runcommandstring = 'canlab_preproc_2012(subjectdir, runwc, imgwc, TR, anat, ''acquisition'', acqorder, ''nointeractive'',''nopart1'',''nodatachecks'',''nomovie'',''nospikeid'',''noorigin'',''nocoreg''); \n';
    jobstr = 'preproc_part2';
    fprintf('Initiating preproc_part2: Warping and smoothing:\n\n')
    for subj = subjs_to_run
        runcommand(subj,EXPT,scriptdir,jobstr,datadir,basedir,runcommandstring,TR,anat,acqorder,runwc,imgwc,dream);
    end
    fprintf('\nPreproc part 2 completed.\n\n')
else
    fprintf('Skipping preproc part 2\n\n')
end


if dream
    fprintf('----------------------------------------------------\n')
    fprintf('|                                                  |\n')
    fprintf('|           Submitting jobs                        |\n')
    fprintf('|                                                  |\n')
    fprintf('----------------------------------------------------\n\n\n')
    
    wd = fullfile(scriptdir,'dream_jobs');
    if exist(wd,'dir'), rmdir(wd,'s'); end
    mkdir(wd)
    addpath(genpath(scriptdir));
       
    % define scheduler
    sched = findResource('scheduler', 'type', 'generic');
    set(sched,'ClusterSize', 1);
    set(sched,'DataLocation',wd);
    set(sched,'ClusterMatlabRoot','/usr/local/matlab/R2011a');
    set(sched,'ParallelSubmitFcn',@parallelSubmitFcn);
    set(sched,'SubmitFcn',@distributedSubmitFcn);
    set(sched,'DestroyJobFcn', @destroyJobFcn);
    set(sched,'DestroyTaskFcn',@destroyTaskFcn);
    
    set(sched,'ClusterOsType','unix');
    set(sched,'HasSharedFileSystem', true);
    nargsout = 0;
    
    scriptlist = filenames(fullfile(scriptdir,'*_publisher.m'));
    if dreamask
        fprintf('Scripts to submit:\n')
        disp(char(scriptlist));        
        while 1;
            resp = input('Proceed with job submission? (y/[n]) ','s');
            if strcmpi(resp,'y')
                break
            else % no
                resp = input('Are you sure? (y/[n]) ','s');
                if strcmpi(resp,'y')
                    fprintf('Exiting\n')
                    return
                end
            end
        end
    end
    
    % submit jobs
    fprintf('%s\tSubmitting jobs to cluster\n',datestr(now,31));
    for s = scriptlist'
        script = regexprep(s{1},'.*/(.*).m','$1');
        j = sched.createJob();
        set(j,'PathDependencies',cellstr(path));
        createTask(j, str2func(script), nargsout, {});
        alltasks = get(j, 'Tasks');
        set(alltasks, 'CaptureCommandWindowOutput', true);
        fprintf('Submitting script: %s\n',script)
        submit(j);
    end
    
    % wait
    % look into using wait, waitForState
    fprintf('WAITING for jobs to stop running (they may run for a while)\n');
    t1 = clock;
    notdone = true;
    fprintf('Elapsed minutes:         ')
    while notdone
        notdone = false;
        for i=1:8, fprintf('\b'); end; fprintf('%8.1f',etime(clock,t1)/60);
        pause(10)
        jobstatefiles = filenames(sprintf('%s/Job*.state.mat',wd),'absolute');
        for s = 1:numel(jobstatefiles)
            jobstate = textread(jobstatefiles{s},'%s');
            if strcmp(jobstate,'running') || strcmp(jobstate,'queued'), notdone = true; break; end
        end
    end
    fprintf('\n%s\tJobs done running\n',datestr(now,31));
    
    % handle output arguments/streams
    fprintf('GATHERING job outputs\n');
    jobdirs = filenames(sprintf('%s/Job*[0-9]',wd),'absolute');
    for i = 1:numel(jobdirs)
        % output stream
        jobout = load(sprintf('%s/Task1.out.mat',jobdirs{i}));
        cw = jobout.commandwindowoutput;
        cw = regexprep(cw,'^> ','');
        cw = regexprep(cw,'\n> ','\n');
        cw = regexprep(cw,'\b[^\n]*\n',' LINES WITH BACKSPACES OMITTED\n');
%         diary(fullfile(wd,sprintf('cmdwnd_%04d.txt',jobin.argsin{2}))),
        disp(cw)
%         diary off
    end
%     % output stream
%     [ignore ignore] = system(sprintf('cat %s/cmdwnd_*txt > %s',wd,fulldiaryname)); %#ok
%     
%     % merge diaries
%     [ignore ignore] = system(sprintf('rm %s/diary*',wd)); %#ok
%     [ignore ignore] = system(sprintf('grep ''^> '' %s | sed ''s|^> ||'' >> %s',fulldiaryname,diaryname)); %#ok
end


fprintf('----------------------------------------------------\n')
fprintf('|                                                  |\n')
fprintf('|           Completed preprocessing                |\n')
fprintf('|                                                  |\n')
fprintf('----------------------------------------------------\n\n\n')

if dream
    fprintf('Submit all the part1 scripts, then reorient all subjets, and then submit\nall part2 scripts.  See wiki for more details\n');
end

end

function runcommand(subj,EXPT,scriptdir,jobstr,datadir,basedir,runcommandstring,TR,anat,acqorder,runwc,imgwc, dream) %#ok
    
    subjectname = EXPT.subjects{subj};
    scriptname = fullfile(scriptdir, ['run_preproc_' subjectname '_' jobstr '.m']);
    
    subjectdir = fullfile(datadir, subjectname);
    
    %scriptname = fullfile(scriptdir, ['run_preproc_' timestr '.m']);
    
    %to make dream compatible
    runcommandstring = strrep(runcommandstring,'subjectdir',['''' subjectdir '''']);
    
    outputdir = fullfile(basedir, 'HTML_output', ['Preprocessing_html_' subjectname]);
    if ~exist(outputdir, 'dir'), mkdir(outputdir), end
    
    p = struct('useNewFigure', false, 'maxHeight', 1500, 'maxWidth', 1200, ...
        'outputDir', outputdir, 'showCode', false);
    
    save(fullfile(scriptdir,'preproc_input.mat'));
    cd(scriptdir)
    
    % - - - - - - - - - - - - - - - - - - - - - - - -
    % BUILD SCRIPT
    % - - - - - - - - - - - - - - - - - - - - - - - -
    % *** use fopen, not diary, to overwrite old scripts if necessary
    %diary(scriptname)
    FID = fopen(scriptname, 'w');
    
    %disp('for subj = subjs_to_run')
    %disp('     subjectname = EXPT.subjects{subj};')
    %disp('     subjectdir = fullfile(datadir, subjectname);');
    
    fprintf(FID, ['%% ' subjectname '\n']);
    fprintf(FID, ['%% Data in: ' subjectdir '\n']);
    fprintf(FID, ['%% Output saved to: ' outputdir '\n']);
    fprintf(FID, ['%% Started on: ' scn_get_datetime '\n']);
    fclose(FID);
    
    eval(['!echo % Created by user: $USER >>' scriptname])
    
    FID = fopen(scriptname, 'a');    
    fprintf(FID, '%% -------------------------------------------------------- \n');
    fprintf(FID, 't1 = clock;\n');
    fprintf(FID, '\n');
    fprintf(FID, 'load preproc_input.mat\n');
%     if dream %must add paths
%         spmpath = '/home/joas2631/spm8';
%         repopath = '/data/projects/wagerlab/RepoYoni';
%         fprintf(FID, ['addpath(genpath(''' repopath '''))\n']);
%         fprintf(FID, ['addpath(genpath(''' spmpath '''))\n']);
%     end
    fprintf(FID, runcommandstring);
    fprintf(FID, '\n');
    %disp('end % for...subject loop');
    %diary off
    
    fprintf(FID, 'elapsed = etime(clock, t1) ./ 60;\n');
    fprintf(FID, 'hours = floor(elapsed ./ 60);\n');
    fprintf(FID, 'min = rem(elapsed, 60);\n');
    fprintf(FID, 'fprintf(''Completed in %%3.0f hours, %%3.0f min\\n'', hours, min);\n');
    
    fclose(FID);
    
    % - - - - - - - - - - - - - - - - - - - - - - - -
    % RUN IT
    % - - - - - - - - - - - - - - - - - - - - - - - -
    fprintf('Preprocessing subject %s\nOutput saved in:%s\n', subjectname, outputdir);
   
    if dream %on dream, write a publisher script to be submitted to cluster
        FID = fopen([scriptname(1:end-2) '_publisher.m'], 'w');
        structstr = ['p = struct(''useNewFigure'', false, ''maxHeight'', 1500, ''maxWidth'', 1200, ''outputDir'',''' outputdir ''', ''showCode'', false);\n']
        fprintf(FID,structstr);
        fprintf(FID, 'publish(''%s'',p)', scriptname);
        fclose(FID);
    else %otherwise, publish now
        publish(scriptname, p);
    end
    
    cd(datadir)
end
