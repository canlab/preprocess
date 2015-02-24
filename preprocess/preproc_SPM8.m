function [PP] = preproc_SPM8(varargin)

% This function is for automated preprocessing of neuroimaging data using
% SPM8, Tor Wager's CANLAB "Core Tools" software (for Matlab), and
% INRIAlign (which should be packaged with this function, do not download 
% existing version from the internet, as they may not work correctly with 
% SPM8). Both SPM8 and the Core Tools code must be in Matlab's path for
% this function to work. It also requires versions of Matlab that support
% empty output variables (use of ~ as an output variable--if you encounter
% errors with code that includes a ~, you can use find/replace in the
% problematic functions to replace ~ with ans and everything should work).
% 
% If the script has been run previously, only data in newly added folders
% will be preprocessed.
% 
% This function does not require all EPI runs to have the same number of
% volumes, nor does it require all subjects to have the same number of
% runs, nor does it have any requirements in terms of filename conventions.
% If it is run once on a dataset, new data can be added to the study and
% the function will "intelligently" only preprocess the newly added data. 
% 
% If preprocessing is interrupted due to an error or system failure,
% rerunning the scripts with switches to disable already-completed
% processing steps allows the user to resume preprocessing on the same set
% of subjects. The script behaves this was automatically, so if the user
% wants to override this behavior they must run the following code from
% their ./study_directory/matlab_code directory: 
% load('PP.mat'); PP = rmfield(PP,'inprogress'); save('PP.mat','PP');
% 
% 
% 
% Carries out the following preprocessing steps in this order (steps named
% in square brackets are not performed by default, all default steps can
% be skipped with the appropriate switches; see below in this help):
% 
% [File type conversion]: Converts all EPI data to 32-bit float.
% 
% [Strip artifacts]: For use on data corrupted by motion-related "striping"
% artifacts. Will replace affected scans from the raw data with the average
% of the unaffected scans temporally adjacent to the affected scans.
% 
% Quality Control: Runs scn_session_spike_id.m on the raw data, to identify
% artifactual scans and produce images for visual checking of data quality
% (see the 'qc_images' folder created for each subject). Also produces a
% mean EPI for each scan to visually check quality of the EPI data.
% 
% Slice Timing: Standard SPM8 slice timing correction. The middle slice is
% used as the reference. Default behavior assumes interleaved slice
% acquisition sequence, but this can be overidden with the 'notinterleaved'
% switch. Creates files with the prefix 'a'.
% 
% Realignment: Motion realignment using INRIAlign. If you are not using 4d
% .img/.hdr files, please manually check that your realignment parameters
% were written to file, as I have not tested INRIAlign's functionality in
% SPM8 with other file types.
% 
% [Check Realignment]: Prints realignment parameters (maximum deviation
% from the median position during each run) and allows the user to exclude
% subjects from further preprocessing based on excessive movement.
% 
% [Pre-coregistration]: Coregisters EPI runs to each other and all EPI runs
% to the structural, prior to manual reorientation.
% 
% Manual Reorientation: Requires the user to manually set the origin and
% orientation of the structural image and the mean EPI image (thus EPI runs
% MUST be in register with each other) to match their respective ICBM
% templates. Matching need not be perfect, this is done to give better
% starting estimates for coregistration and normalization.
% 
% Coregistration: Coregisters the structural image to the ICBM template,
% and coregisters all EPI runs to the structural.
% 
% Check Registration: Requires the user to visually verify that EPI images
% are in register with the structural for each subject.
% 
% Segmentation: Segments the structural image into the three tissue
% compartments. MUST be completed prior to normalization, as the
% normalization parameters are taken from the segmentation results.
% 
% Normalization: Normalizes structural image to the ICBM template and
% applies the normalization to all EPI images. Creates files with the
% prefix 'wa'.
% 
% Check Normalization: Requires the user to visually verify that the
% normalization was successful. 
% 
% Smoothing: Smooths images with an 8mm FWHM kernel. At present this cannot
% be changed, but this functionality should be added soon. Creates files
% with the prefix 'swa'.
% 
% [Normalize image values]: Divides each voxel by the mean of that voxel's
% timeseries and multiplies by 100. This is done to make first-level HRF
% response estimates (not betas) scale to % signal change after first-level
% modelling. Creates files with the prefix 'nswa'.
% 
% 
% 
% Usage: PP = preproc
%
% This usage runs the script in the pwd, assuming each existing directory
% starting with the letter 'S' contains data for a separate subject (the 
% 'S' prefix can be changed to something else by the user), with
% structurals in the subject's root directory and EPIs in a folder called 
% "EPIs".
% 
% That is, your directory structure MUST have the following format:
% ./STUDY_BASEDIRECTORY
% ./STUDY_BASEDIRECTORY/SUBJECT1_DIRECTORY
% ./STUDY_BASEDIRECTORY/SUBJECT1_DIRECTORY/structural_image_filename
% ./STUDY_BASEDIRECTORY/SUBJECT1_DIRECTORY/EPIs
% ./STUDY_BASEDIRECTORY/SUBJECT1_DIRECTORY/EPIs/functional_image_filenames
% ./STUDY_BASEDIRECTORY/SUBJECT2_DIRECTORY
% and so on...
%
% Usage: PP = preproc(PP)
%
% This allows specification of various fields to change function behavior:
% Allowed fields:
%     basedir:     <string> Full path to the experiment's base directory.
%                  Allows the user to call the function from any directory.
%     TR:          <scalar> Your study's TR in seconds. You will prompted
%                  for this info if it is not provided.
%     forcedo:     <vector> A list of subject numbers (corresponding to the
%                  elements of PP.dirs in ./matlab_code/PP.mat). This
%                  forces the script to redo preprocessing on these
%                  subjects, instead of looking only for new subjects to
%                  preprocess.
%     MADs:        <scalar> Changes the threshold, in Mean Absolute
%                  Deviations, during the QC procedure for identifying
%                  "bad" scans (default = 10).
%     excl:        <vector> A list of subject numbers (corresponding to the
%                  order of directories in PP.basedir or pwd) to *not* run
%                  preprocessing on.
%     prefix:      <string> The string that all subject directories start
%                  with. Directory names that do not begin with this string
%                  will be assumed not to contain imaging data for
%                  preprocessing (default = 'S').
%     segment_job: A cell structure in SPM8 matlabbatch format for
%                  segmentation (segment_job{1}.spm.spatial.preproc).
%                  Entered values will override the scripts default
%                  segmentation.
%     smooth_job:  A cell structure in SPM8 matlabbatch format for
%                  smoothing (segment_job{1}.spm.spatial.smooth). Entered
%                  values will override the scripts default smoothing.
%                  Mostly useful if you set
%                  PP.segment_job{1}.spm.spatial.smooth.fwhm = [X X X];
%                  where X is your kernel size in MM. Defaults to [8 8 8].
% 
% Eventually I intend to allow overwriting of default inputs to all of the
% SPM8 jobs called by this function. If you have a specific request please
% let me know and I'll try to implement it ASAP.
% 
% 
% Usage: PP = preproc(PP,switch1,switch2,...)
% 
% Several string switches can be used to change script behavior. They may
% be specified in any order. See above for more details on each
% preprocessing step. Allowed switches:
% 
% 'nifti'          Causes the script to look for and work on '.nii' files
%                  instead of '.img' files. PLEASE NOTE: I have not tested
%                  the motion realignment with nifti images. If you're
%                  using nifti files please ensure that after realignment
%                  the values were actually written to the images (ask me
%                  if you need to know how to do this).
% 
% 'skipdefaults'   Use if you have already run the script during this
%                  instance of matlab. Not required, but will save time on
%                  slower machines by skipping calls to
%                  spm('defaults,'FMRI') and spm_jobman('initcfg').
% 
% 'doconvert'      Enables conversion of images to 32-bit float.
% 
% 'stripartifacts' Removes motion-related 'striping' artifacts from the raw
%                  data.
% 
% 'noQC'           Disables quality control measures.
% 
% 'notinterleaved' Use this if your data was not acquired with an
%                  interleaved acquisition sequence.
%
% 'notPhilipsinterleaved' Use this if your data was interleaved, but not on
%                         a Philips scanner. This is for acquisition
%                         sequences that go [1:2:end 2:2:end]. The default
%                         behavior is [1:round(sqrt(num_slices)):end
%                         2:round(sqrt(num_slices)):end ... etc.
% 
% 'noslicetiming'  Disables slice timing correction.
% 
% 'norealign'      Disables motion realignment.
% 
% 'checkrealign'   Prints realignment parameters, even if realignment was
%                  not performed on this run of the function.
% 
% 'realigndiff'    Reports the maximum derivative of each realignment
%                  parameter instead of the default behavior of reporting
%                  the maximum absolute deviation from the median position
%                  of the scan.
% 
% 'doprecoreg'     Enables pre-coregistration.
%  
% 'noreorient'     Disables manual reorientation.
% 
% 'nocoreg'        Disables coregistration.
% 
% 'checkreg'       Use only with 'nocoreg'. This will allow the user to
%                  check registration between the structural and EPI images
%                  even if coregistration was not performed.
% 
% 'nosegment'      Disables segmentation.
% 
% 'nonorm'         Disables normalization.
% 
% 'checknorm'      Use only with 'nonorm'. This will allow the user to
%                  check normalization of the structural and EPI images to
%                  the ICBM template.
% 
% 'nosmooth'       Disables smoothing.
% 
% 'dowritemean'    Writes a "mean_functional.img" from the mean of all
%                  subjects' normalized structural images into the
%                  matlab_code directory. 
% 
% 'donormimages'   Normalizes the values in each voxel with respect to the
%                  mean in the time course (see description above).
% 

% Copyright Jared Van Snellenberg 2011. Permission granted for all
% non-commercial uses.



if isempty(varargin) || ~isstruct(varargin{1})
    PP.basedir = pwd;
else
    PP = varargin{1};
    if ~isfield(PP,'basedir')
        PP.basedir = pwd;
    end
end

segment_job{1}.spm.spatial.preproc.output.GM = [0 1 1];
segment_job{1}.spm.spatial.preproc.output.WM = [0 1 1];
segment_job{1}.spm.spatial.preproc.output.CSF = [0 0 0];
segment_job{1}.spm.spatial.preproc.output.biascor = 1;
segment_job{1}.spm.spatial.preproc.output.cleanup = 0;

segment_job{1}.spm.spatial.preproc.opts.tpm = {which(['tpm' filesep 'grey.nii']),which(['tpm' filesep 'white.nii']),which(['tpm' filesep 'csf.nii'])}';
segment_job{1}.spm.spatial.preproc.opts.ngaus = [2; 2; 2; 4];
segment_job{1}.spm.spatial.preproc.opts.regtype = 'mni';
segment_job{1}.spm.spatial.preproc.opts.warpreg = 1;
segment_job{1}.spm.spatial.preproc.opts.warpco = 20;
segment_job{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
segment_job{1}.spm.spatial.preproc.opts.biasfwhm = 60;
segment_job{1}.spm.spatial.preproc.opts.samp = 2;
segment_job{1}.spm.spatial.preproc.opts.mask = {''};

if isfield(PP,'segment_job')
    v = PP.segment_job{1}.spm.spatial.preproc;
    if isfield(v,'output')
        c = struct2cell(v.output);
        f = fieldnames(v.output);
        for i = 1:length(f)
            segment_job{1}.spm.spatial.preproc.output.(f{i}) = c{i};
        end
    end
    if isfield(v,'opts')
        c = struct2cell(v.opts);
        f = fieldnames(v.opts);
        for i = 1:length(f)
            segment_job{1}.spm.spatial.preproc.opts.(f{i}) = c{i};
        end
    end
end

smooth_job{1}.spm.spatial.smooth.fwhm = [8 8 8];
smooth_job{1}.spm.spatial.smooth.dtype = 0;
smooth_job{1}.spm.spatial.smooth.im = 0;
smooth_job{1}.spm.spatial.smooth.prefix = 's';
            
if isfield(PP,'smooth_job')
    v = PP.smooth_job{1}.spm.spatial.smooth;
    c = struct2cell(v);
    f = fieldnames(v);
    for i = 1:length(f)
        smooth_job{1}.spm.spatial.smooth.(f{i}) = c{i};
    end
end


doconvert = 0;
stripartifacts = 0;
doslicetiming = 1;
dorealign = 1;
checkrealign = 0;
skipdefaults = 0;
doQC = 1;
doprecoreg = 0;
docoreg = 1;
doreorient = 1;
checkreg = 1;
forcereg = 0;
checknorm = 1;
dosegment = 1;
dosmooth = 1;
donorm = 1;
forcenorm = 0;
writemean = 0;
ftype = '.img';
realigndiff = 0;
normimages = 0;
interleaved = 1;

MADs = 10;
if isfield(PP,'MADs')
    MADs = PP.MADs;
end

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch(varargin{k})
            case 'doconvert'
                doconvert = 1;
            case 'stripartifacts'
                stripartifacts = 1;
            case 'noslicetiming'
                doslicetiming = 0;
            case 'notinterleaved'
                interleaved = 0;
            case 'notPhilipsinterleaved'
                interleaved = 2;
            case 'norealign'
                dorealign = 0;
            case 'checkrealign'
                checkrealign = 1;
            case 'skipdefaults'
                skipdefaults = 1;
            case 'noQC'
                doQC = 0;
            case 'doprecoreg'
                doprecoreg = 1;
            case 'nocoreg'
                docoreg = 0;
                if ~forcereg
                    checkreg = 0;
                end
            case 'noreorient'
                doreorient = 0;
            case 'checkreg'
                checkreg = 1;
                forcereg = 1;
            case 'nosegment'
                dosegment = 0;
            case 'nosmooth'
                dosmooth = 0;
            case 'nonorm'
                donorm = 0;
                if ~forcenorm
                    checknorm = 0;
                end
            case 'checknorm'
                forcenorm = 1;
                checknorm = 1;
            case 'dowritemean'
                writemean = 1;
            case 'nifti'
                ftype = '.nii';
            case 'realigndiff'
                realigndiff = 1;
            case 'donormimages'
                normimages = 1;
        end
    end
end

if ~skipdefaults
    spm('defaults','FMRI');
    spm_jobman('initcfg');
end
PWD = pwd;

if ~isfield(PP,'basedir')
    PP.basedir = pwd;
end

if isfield(PP,'forcedo')
    forcedo = PP.forcedo;
end

if isfield(PP,'prefix')
    prefix = PP.prefix;
end

if ~isdir([PP.basedir filesep 'matlab_code'])
    mkdir([PP.basedir filesep 'matlab_code']);
end
cd([PP.basedir filesep 'matlab_code']);
if isfield(PP,'excl');
    e = PP.excl;
end
if exist('PP.mat')
    disp('Loading existing PP structure from saved .mat file')
    load PP.mat PP
end
if exist('e')==1
    PP.excl = e;
end

if exist('prefix') == 1
    PP.prefix = prefix;
elseif ~isfield(PP,'prefix')
    PP.prefix = 'S';
end

if exist('forcedo') ~= 1
    todo = [];
    if ~isfield(PP,'dirs')
        dirs = dir([PP.basedir filesep PP.prefix '*']);
        r = [];
        for i = 1:length(dirs)
            r(end+1) = ~dirs(i).isdir;
        end
        dirs(logical(r)) = [];
        todo = 1:length(dirs);
    elseif length(dir([PP.basedir filesep PP.prefix '*'])) > length(PP.dirs)
        if ~isfield(PP,'excl') || length(dir([PP.basedir filesep PP.prefix '*'])) > length(PP.dirs) + length(PP.excl)
            oldirs = PP.dirs;
            dirs = dir([PP.basedir filesep PP.prefix '*']);
            for i = 1:length(dirs)
                new = 1;
                for j = 1:length(oldirs)
                    if strcmp(dirs(i).name,oldirs(j).name)
                        new = 0;
                        break
                    end
                end
                if new
                    todo(end+1) = i;
                end
            end
        else
            disp('Nothing to do!');
            cd(PWD);
            return
        end
    else
        disp('Nothing to do!');
        cd(PWD);
        return
    end
else
    if ischar(forcedo) && strcmp(forcedo,'all')
        forcedo = 1:length(PP.dirs);
    end
    todo = forcedo;
    todo = sort(todo);
    images = PP.images(todo,:);
    dirs = PP.dirs(todo);
    structurals = PP.structurals(todo);
    mfiles = PP.mfiles(todo,:);
end

if isfield(PP,'excl') && exist('forcedo') ~= 1
    for i = 1:length(PP.excl)
        todo(todo==PP.excl(i)) = [];
    end
else
    PP.excl = [];
end

if exist('forcedo') ~= 1
    dirs = dirs(todo);
end
    
if isfield(PP,'inprogress')
    disp('An inprogress analysis was terminated early, loading filename variables')
    ip = PP.inprogress;
    if isfield(ip,'structurals')
        structurals = ip.structurals;
    end
    if isfield(ip,'images')
        images = ip.images;
    end
    if isfield(ip,'dirs')
        dirs = ip.dirs;
    end
    if isfield(ip,'mfiles')
        mfiles = ip.mfiles;
    end
    if isfield(ip,'afiles')
        afiles = ip.afiles;
    end
    if isfield(ip,'wstructurals')
        wstructurals = ip.wstructurals;
    end
    if isfield(ip,'wafiles')
        wafiles = ip.wafiles;
    end
    if isfield(ip,'mean_func')
        mean_func = ip.mean_func;
    end
    if isfield(ip,'swafiles')
        swafiles = ip.swafiles;
    end
end

if ~isfield(PP,'TR')
    PP.TR = input(['Please enter your TR in seconds:\n']);
end
TR = PP.TR;

if exist('images')~=1
    for i = 1:length(dirs)
        t{i} = filenames([PP.basedir filesep dirs(i).name filesep '*' ftype], 'char', 'absolute');
        structurals{i} = strrep(t{i}(1,:),' ','');
        for j = 1:size(t{i},1)-1
            images{i,j} = deblank(t{i}(j+1,:));
        end
    end
end

disp('Running preprocessing steps on the following subjects:');
for i = 1:length(dirs)
    disp(dirs(i).name);
end

PP.inprogress.dirs = dirs;
PP.inprogress.images = images;
PP.inprogress.structurals = structurals;
save PP PP

if doconvert
    jochen_convert(images);
end

if stripartifacts
    bad = strip_motion_artifacts(images);
    for i = 1:size(bad,1)
        badscans = bad(i,:);
        save([PP.basedir filesep dirs(i).name filesep 'bad_scans.mat'],'badscans');
    end
end

if doQC
    for i = 1:length(dirs)
        cd(['..' filesep dirs(i).name])
        scn_session_spike_id(stripemptycells(images(i,:)),[],MADs);
        
        for j = 1:size(images,2)
            if ~isempty(images{i,j})
                V = spm_vol(images(i,j));
                reslice_mean{1}.spm.spatial.realign.write.roptions.which = [0 1];
                reslice_mean{1}.spm.spatial.realign.write.roptions.interp = 4;
                reslice_mean{1}.spm.spatial.realign.write.roptions.wrap= [0 0 0];
                reslice_mean{1}.spm.spatial.realign.write.roptions.mask = 1;
                reslice_mean{1}.spm.spatial.realign.write.roptions.prefix = 'r';
                a = spm5_image_list(length(V{1}), images(i,j));
                reslice_mean{1}.spm.spatial.realign.write.data = a{1};
                spm_jobman('run',{reslice_mean});
                f = strfind(images{i,j},filesep);
                spm_image('init',[images{i,j}(1:f(end)) 'mean' images{i,j}(f(end)+1:end)]);
                saveas(gcf,[PP.basedir filesep dirs(i).name filesep 'qc_images' filesep 'mean_EPI_run' num2str(j) '.png']);
            end
        end
    end
    cd(['..' filesep 'matlab_code'])
    close all
end




for i = 1:size(images,1)
    for j = 1:size(images,2)
        if ~isempty(images{i,j})
            func_files = images(i,j);
            V = spm_vol(func_files);
            num_slices = V{1}(1).dim(3);
            if doslicetiming
                slice_timing_job{1}.spm.temporal.st.scans = spm5_image_list(length(V{1}), func_files);
                slice_timing_job{1}.spm.temporal.st.nslices = num_slices;
                slice_timing_job{1}.spm.temporal.st.tr = TR;
                slice_timing_job{1}.spm.temporal.st.ta = TR - (TR / num_slices);
                if interleaved==1
                    slice_timing_job{1}.spm.temporal.st.so = [];
                    for k = 1:round(sqrt(num_slices))
                        slice_timing_job{1}.spm.temporal.st.so = [slice_timing_job{1}.spm.temporal.st.so k:round(sqrt(num_slices)):num_slices];
                    end
                elseif interleaved==2
                    slice_timing_job{1}.spm.temporal.st.so = [1:2:num_slices 2:2:num_slices]; %interleaved acquisition
                else
                    slice_timing_job{1}.spm.temporal.st.so = 1:num_slices; %not interleaved
                end
                slice_timing_job{1}.spm.temporal.st.refslice = round(num_slices/2);
                slice_timing_job{1}.spm.temporal.st.prefix = 'a';
                
                spm_jobman('run', {slice_timing_job});
            end
            
            f = strfind(func_files{1},filesep);
            afiles{i,j} = [func_files{1}(1:f(end)) 'a' func_files{1}(f(end)+1:end)];
            
            PP.inprogress.afiles = afiles;
            save PP PP
            
            if dorealign
                inria_realign(afiles{i,j}, struct('quality',1,'fwhm',8));
            
                reslice_mean{1}.spm.spatial.realign.write.roptions.which = [0 1];
                reslice_mean{1}.spm.spatial.realign.write.roptions.interp = 4;
                reslice_mean{1}.spm.spatial.realign.write.roptions.wrap= [0 0 0];
                reslice_mean{1}.spm.spatial.realign.write.roptions.mask = 1;
                reslice_mean{1}.spm.spatial.realign.write.roptions.prefix = 'r';
                a = spm5_image_list(length(V{1}), afiles(i,j));
                reslice_mean{1}.spm.spatial.realign.write.data = a{1};
                
                spm_jobman('run',{reslice_mean});
            end
            
            f = strfind(afiles{i,j},filesep);
            mfiles{i,j} = [afiles{i,j}(1:f(end)) 'mean' afiles{i,j}(f(end)+1:end)];
            
            PP.inprogress.mfiles = mfiles;
            save PP PP
        end
    end
end


if checkrealign
    clear d
    clear c
    
    for i = 1:size(images,1)
        d{i} = filenames(['..' filesep dirs(i).name filesep '*.txt'], 'absolute');
    end
    
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if ~isempty(images{i,j})
                rp{i,j} = load(d{i}{j});
            end
        end
    end
    
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if ~isempty(images{i,j})
                if realigndiff
                    for k = 1:6
                        dRP{k} = diff(rp{i,j}(:,k));
                    end
                else
                    for k = 1:6
                        dRP{k} = rp{i,j}(:,k) - median(rp{i,j}(:,k));
                    end
                end
                fprintf(1,['Index ' num2str(i) ': ' dirs(i).name ' Run: ' num2str(j) ' Max RPs:\t'...
                    num2str(max(abs(dRP{1}))) '\t' ...
                    num2str(max(abs(dRP{2}))) '\t' ...
                    num2str(max(abs(dRP{3}))) '\t' ...
                    num2str(max(abs(dRP{4}))*180/pi) '\t' ...
                    num2str(max(abs(dRP{5}))*180/pi) '\t' ...
                    num2str(max(abs(dRP{6}))*180/pi) '\t' ...
                    '\n']);
            end
        end
    end
    
    excl = [];
    resp = 1;
    beep
    while ~isempty(resp)
        resp = input('Check the printed realignment parameters. Enter an index number to exclude that subject!\n');
        if ~isempty(resp)
            excl(end+1) = resp;
        end
    end
    images(excl,:) = [];
    afiles(excl,:) = [];
    mfiles(excl,:) = [];
    structurals(excl) = [];
    dirs(excl) = [];
    if isfield(PP,'excl')
        PP.excl = [PP.excl excl];
    else
        PP.excl = excl;
    end
end

try PP.inprogress.images = images; catch, end
try PP.inprogress.afiles = afiles; catch, end
try PP.inprogress.mfiles = mfiles; catch, end
try PP.inprogress.structurals = structurals; catch, end
try PP.inprogress.dirs = dirs; catch, end
save PP PP


if doprecoreg && docoreg
    precoreg_job{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    precoreg_job{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    precoreg_job{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.0200    0.0200    0.0200    0.0010    0.0010    0.0010    0.0100    0.0100    0.0100    0.0010    0.0010    0.0010];
    precoreg_job{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            precoreg_job{i,j} = precoreg_job{1};
            if j ~= size(images,2)
                precoreg_job{i,j}.spm.spatial.coreg.estimate.ref = {[mfiles{i,size(images,2)-j} ',1']};
                precoreg_job{i,j}.spm.spatial.coreg.estimate.source = {[mfiles{i,size(images,2)-j+1} ',1']};
                list = cell({});
                for k = 1:j
                    if ~isempty(images{i,size(images,2)-k+1})
                        V = spm_vol(afiles(i,size(images,2)-k+1));
                        a = spm5_image_list(length(V{1}), afiles(i,size(images,2)-k+1));
                        list(end+1:end+length(V{1})) = a{1};
                    end
                end
                if isempty(list)
                    precoreg_job{i,j} = [];
                    continue
                end
                precoreg_job{i,j}.spm.spatial.coreg.estimate.other = list';
            else
                precoreg_job{i,j}.spm.spatial.coreg.estimate.ref = {[structurals{i} ',1']};
                precoreg_job{i,j}.spm.spatial.coreg.estimate.source = {[mfiles{i,1} ',1']};
                list = cell({});
                for k = 1:j
                    if ~isempty(images{i,size(images,2)-k+1})
                        V = spm_vol(afiles(i,size(images,2)-k+1));
                        a = spm5_image_list(length(V{1}), afiles(i,size(images,2)-k+1));
                        list(end+1:end+length(V{1})) = a{1};
                    end
                end
                precoreg_job{i,j}.spm.spatial.coreg.estimate.other = list';
            end
        end
    end
    
    
    spm_jobman('run',{stripemptycells(precoreg_job(:)')});
    
    
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if ~isempty(images{i,j})
                V = spm_vol(afiles(i,j));
                reslice_mean{1}.spm.spatial.realign.write.roptions.which = [0 1];
                reslice_mean{1}.spm.spatial.realign.write.roptions.interp = 4;
                reslice_mean{1}.spm.spatial.realign.write.roptions.wrap= [0 0 0];
                reslice_mean{1}.spm.spatial.realign.write.roptions.mask = 1;
                reslice_mean{1}.spm.spatial.realign.write.roptions.prefix = 'r';
                a = spm5_image_list(length(V{1}), afiles(i,j));
                reslice_mean{1}.spm.spatial.realign.write.data = a{1};
                
                spm_jobman('run',{reslice_mean});
            end
        end
    end
end


toreorient = 1:size(images,1);

while ~isempty(toreorient)
    
    if doreorient
        for i = toreorient
            spm_image('init',structurals{i});
            input(['Reorient this T1 and apply.\nImage should be manually reoriented to match which(''single_subj_T1.nii'')\n Currently viewing ' dirs(i).name]);
            if ~isempty(images{i,1})
                spm_image('init',mfiles{i,1});
            else
                for j = 1:size(images,2)
                    if ~isempty(images{i,j})
                        spm_image('init',mfiles{i,j});
                        break
                    end
                end
            end
            input(['Reorient this mean image and apply to all functionals.\nImage should be manually reoriented to match which(''EPI.nii'')\n Currently viewing ' dirs(i).name]);
        end
    end
    
    if docoreg
        coreg_defaults{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        coreg_defaults{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        coreg_defaults{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.0200    0.0200    0.0200    0.0010    0.0010    0.0010    0.0100    0.0100    0.0100    0.0010    0.0010    0.0010];
        coreg_defaults{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        
        for i = toreorient
            for j = 1:size(images,2) - 1
                if ~isempty(images{i,j})
                    coreg_job{i,j} = coreg_defaults{1};
                    coreg_job{i,j}.spm.spatial.coreg.estimate.ref = {[mfiles{i,size(images,2)-j} ',1']};
                    coreg_job{i,j}.spm.spatial.coreg.estimate.source = {[mfiles{i,size(images,2)-j+1} ',1']};
                    list = cell({});
                    for k = 1:j
                        if ~isempty(images{i,size(images,2)-k+1})
                            V = spm_vol(afiles(i,size(images,2)-k+1));
                            a = spm5_image_list(length(V{1}), afiles(i,size(images,2)-k+1));
                            list(end+1:end+length(V{1})) = a{1};
                        end
                    end
                    if isempty(list)
                        coreg_job{i,j} = [];
                        continue
                    end
                    coreg_job{i,j}.spm.spatial.coreg.estimate.other = list';
                end
            end
            if i == toreorient(1) 
                if size(images,2) ~= 1
                    coreg_job{i,end+1} = coreg_defaults{1};
                else
                    coreg_job{i,1} = coreg_defaults{1};
                end
                coreg_job{i,end}.spm.spatial.coreg.estimate.ref = {which('T1.nii')};
                coreg_job{i,end}.spm.spatial.coreg.estimate.source = {[structurals{i} ',1']};
                coreg_job{i,end}.spm.spatial.coreg.estimate.other = cell({});
                coreg_job{i,end}.spm.spatial.coreg.estimate = rmfield(coreg_job{i,end}.spm.spatial.coreg.estimate,'other');
                
                coreg_job{i,end+1} = coreg_defaults{1};
            else
                coreg_job{i,end-1} = coreg_defaults{1};
                coreg_job{i,end-1}.spm.spatial.coreg.estimate.ref = {which('T1.nii')};
                coreg_job{i,end-1}.spm.spatial.coreg.estimate.source = {[structurals{i} ',1']};
                coreg_job{i,end-1}.spm.spatial.coreg.estimate.other = cell({});
                coreg_job{i,end-1}.spm.spatial.coreg.estimate = rmfield(coreg_job{i,end-1}.spm.spatial.coreg.estimate,'other');
                coreg_job{i,end} = coreg_defaults{1};
            end
            coreg_job{i,end}.spm.spatial.coreg.estimate.ref = {[structurals{i} ',1']};
            coreg_job{i,end}.spm.spatial.coreg.estimate.source = {[mfiles{i,1} ',1']};
            list = cell({});
            for k = 1:size(images,2)
                if ~isempty(images{i,size(images,2)-k+1})
                    V = spm_vol(afiles(i,size(images,2)-k+1));
                    a = spm5_image_list(length(V{1}), afiles(i,size(images,2)-k+1));
                    list(end+1:end+length(V{1})) = a{1};
                end
            end
            coreg_job{i,end}.spm.spatial.coreg.estimate.other = list';
        end
        
        spm_jobman('run',{stripemptycells(coreg_job(:)')});
        
        for i = toreorient
            for j = 1:size(images,2)
                if ~isempty(images{i,j})
                    V = spm_vol(afiles(i,j));
                    reslice_mean{1}.spm.spatial.realign.write.roptions.which = [0 1];
                    reslice_mean{1}.spm.spatial.realign.write.roptions.interp = 4;
                    reslice_mean{1}.spm.spatial.realign.write.roptions.wrap= [0 0 0];
                    reslice_mean{1}.spm.spatial.realign.write.roptions.mask = 1;
                    reslice_mean{1}.spm.spatial.realign.write.roptions.prefix = 'r';
                    a = spm5_image_list(length(V{1}), afiles(i,j));
                    reslice_mean{1}.spm.spatial.realign.write.data = a{1};
                    
                    spm_jobman('run',{reslice_mean});
                end
            end
        end
    end
    
    if checkreg
        doagain = [];
        for i = toreorient
            l = length(stripemptycells(mfiles(i,:)));
            for j = 1:ceil(l/5)
                spm_check_registration(strvcat(structurals{i},mfiles{i,(j-1)*5+1:min([(j-1)*5+5 l])}));
                didit = 0;
                while ~didit
                    try
                        doagain(end+1) = input(['Checking registration for ' dirs(i).name ' Runs: ' num2str((j-1)*5+1) ' to ' num2str(min([(j-1)*5+5 l])) ', enter 0 to go to the next subject\n       OR ENTER ' num2str(i) ' TO REDO THIS SUBJECT:\n']);
                        didit = 1;
                    catch
                        didit = 0;
                    end
                end
            end
        end
        toreorient = unique(doagain);
        toreorient(toreorient==0) = [];
        clear coreg_job
    else
        toreorient = [];
    end
end


if dosegment   
    for i = 1:length(structurals)
        segment_job{1}.spm.spatial.preproc.data{i,1} = [structurals{i} ',1'];
    end
    spm_jobman('run',{segment_job});
end

if donorm
    write_norm{1}.spm.spatial.normalise.write.roptions.preserve = 0;
    write_norm{1}.spm.spatial.normalise.write.roptions.bb = [-78 -120 -80; 78 95 140];
    write_norm{1}.spm.spatial.normalise.write.roptions.vox = [3 3 3];
    write_norm{1}.spm.spatial.normalise.write.roptions.interp = 1;
    write_norm{1}.spm.spatial.normalise.write.roptions.wrap = [0 1 0];
    write_norm{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
    write_norm{2} = write_norm{1};
    write_norm{2}.spm.spatial.normalise.write.roptions.vox = [1 1 1];
    
    for i = 1:size(images,1)
        write_norm{1}.spm.spatial.normalise.write.subj.matname = {strrep(structurals{i},ftype,'_seg_sn.mat')};
        write_norm{2}.spm.spatial.normalise.write.subj.matname = write_norm{1}.spm.spatial.normalise.write.subj.matname;
        
        
        list = cell({});
        for j = 1:size(images,2)
            if ~isempty(images{i,j})
                V = spm_vol(afiles(i,j));
                a = spm5_image_list(length(V{1}), afiles(i,j));
                list(end+1:end+length(V{1})) = a{1};
            end
        end
        write_norm{1}.spm.spatial.normalise.write.subj.resample = list';
        
        write_norm{2}.spm.spatial.normalise.write.subj.resample = structurals(i);
        
        spm_jobman('run',{write_norm});
    end
end


for i = 1:size(images,1)
    for j = 1:size(images,2)
        if ~isempty(images{i,j})
            A = findstr(filesep,images{i,j});
            wafiles{i,j} = [images{i,j}(1:A(end)) 'wa' images{i,j}(A(end)+1:end)];
        end
    end
    A = findstr(filesep,structurals{i});
    wstructurals{i} = [structurals{i}(1:A(end)) 'w' structurals{i}(A(end)+1:end)];
end

PP.inprogress.wafiles = wafiles;
PP.inprogress.wstructurals = wstructurals;
save PP PP

for i = 1:size(images,1)
    list = cell('');
    for j = 1:size(images,2)
        if ~isempty(images{i,j})
            list(j,1) = wafiles(i,j);
        end
    end
    if donorm
        mean_image(list,['..' filesep dirs(i).name filesep 'EPIs' filesep 'mean_warped_functional' ftype]);
    end
    mean_func{i} = filenames(['..' filesep dirs(i).name filesep 'mean_warped_functional' ftype], 'char', 'absolute');
end

if checknorm
    for i = 1:size(images,1)
        spm_check_registration(strvcat(which('T1.nii'),wstructurals{i},mean_func{i}));
        input(['Checking registration for ' dirs(i).name]);
    end
end

PP.inprogress.mean_func = mean_func;
save PP PP

if dosmooth
    for i = 1:size(images,1)
        list = cell({});
        for j = 1:size(images,2)
            if ~isempty(images{i,j})
                V = spm_vol(wafiles(i,j));
                a = spm5_image_list(length(V{1}), wafiles(i,j));
                list(end+1:end+length(V{1})) = a{1};
            end
        end
        smooth_job{1}.spm.spatial.smooth.data = list';
        
        spm_jobman('run',{smooth_job});
    end
end


for i = 1:size(images,1)
    for j = 1:size(images,2)
        if ~isempty(images{i,j})
            A = findstr(filesep,images{i,j});
            swafiles{i,j} = [images{i,j}(1:A(end)) 'swa' images{i,j}(A(end)+1:end)];
        end
    end
end

PP.inprogress.swafiles = swafiles;
save PP PP

if writemean
    mean_image(wstructurals,['mean_warped_structural' ftype]);
end

if normimages
    for i = 1:size(swafiles,1)
        for j = 1:size(swafiles,2)
            if ~isempty(swafiles{i,j})
                A = findstr(filesep,images{i,j});
                nswafiles{i,j} = [images{i,j}(1:A(end)) 'nswa' images{i,j}(A(end)+1:end)];
                
                V = spm_vol(swafiles{i,j});
                Y = spm_read_vols(V);
                y = Y;
                for k = 1:size(Y,4)
                    y(:,:,:,k) = 100 * Y(:,:,:,k) ./ mean(Y,4);
                    V(k).fname = nswafiles{i,j};
                    spm_write_vol(V(k),y(:,:,:,k));
                end
            end
        end
    end
end

%structurals images dirs mfiles afiles wstructurals wafiles mean_func fmask smask swafiles


if ~isfield(PP,'images')
    PP.dirs = dirs;
    PP.structurals = structurals;
    PP.images = images;
    PP.mfiles = mfiles;
    PP.afiles = afiles;
    PP.wstructurals = wstructurals;
    PP.wafiles = wafiles;
    PP.mean_func = mean_func;
    PP.swafiles = swafiles;
    if exist('nswafiles') == 1
        PP.nswafiles = nswafiles;
    end
% elseif exist('forcedo') ~= 1 && length(todo)+length(PP.structurals)-length(PP.excl) ~= length(PP.dirs)
%     error('Something''s up with the filename variables, please fix in debug mode')
else %%this part doesn't quite work right yet, but the above error should catch it when it's going to break
    if exist('forcedo') ~= 1
        if isfield(PP,'excl')
            indx = 1:length(dir([PP.basedir filesep PP.prefix '*']))-length(PP.excl);
            for i = 1:length(PP.excl)
                if any(PP.excl > todo)
                    error('This is a total kludge, I haven''t worked out a decent solution yet. If you want this script to work correctly in the future you need to set the ''indx'' variable to a vector of all the OLD subjects in your study, and the ''todo'' variable to a vector of all the new ones, and then run all of the code that occurs after this error text until the main function end')
                end
            end
            todo = todo - length(PP.excl);
        else
            indx = 1:length(dir([PP.basedir filesep PP.prefix '*']));
        end
        indx(todo) = [];
        
        
        PP.dirs(indx) = PP.dirs;
        PP.structurals(indx) = PP.structurals;
        PP.images(indx,:) = PP.images;
        PP.mfiles(indx,:) = PP.mfiles;
        PP.afiles(indx,:) = PP.afiles;
        PP.wstructurals(indx) = PP.wstructurals;
        PP.wafiles(indx,:) = PP.wafiles;
        PP.mean_func(indx) = PP.mean_func;
        PP.swafiles(indx,:) = PP.swafiles;
        if exist('nswafiles') == 1
            PP.nswafiles(indx,:) = PP.nswafiles;
        end
        
        PP.dirs(todo) = dirs;
        if size(images,2) > size(PP.images,2)
            PP.images(:,size(images,2)) = cell(1);
            PP.mfiles(:,size(images,2)) = cell(1);
            PP.afiles(:,size(images,2)) = cell(1);
            PP.wafiles(:,size(images,2)) = cell(1);
            PP.swafiles(:,size(images,2)) = cell(1);
            PP.nswafiles(:,size(images,2)) = cell(1);
        elseif size(PP.images,2) > size(images,2)
            images(:,size(PP.images,2)) = cell(1);
            mfiles(:,size(PP.images,2)) = cell(1);
            afiles(:,size(PP.images,2)) = cell(1);
            wafiles(:,size(PP.images,2)) = cell(1);
            swafiles(:,size(PP.images,2)) = cell(1);
            nswafiles(:,size(PP.images,2)) = cell(1);
        end
        PP.structurals(todo) = structurals;
        PP.images(todo,:) = images;
        PP.mfiles(todo,:) = mfiles;
        PP.afiles(todo,:) = afiles;
        PP.wstructurals(todo) = wstructurals;
        PP.wafiles(todo,:) = wafiles;
        PP.mean_func(todo) = mean_func;
        PP.swafiles(todo,:) = swafiles;
        if exist('nswafiles') == 1
            PP.nswafiles(todo,:) = nswafiles;
        end
        
    end
end

if isfield(PP,'inprogress')
    PP = rmfield(PP,'inprogress');
end
save PP PP

cd(PWD)
end


function [] = jochen_convert(images)
for i = 1:size(images,1)
    for j = 1:size(images,2)
        if ~isempty(images{i,j})
            % get volume descriptor (make sure to call spm('defaults', 'FMRI') first!)
            v = spm_vol(images{i,j});
            
            % read data (into 4D variable, double precision internally)
            y = spm_read_vols(v);
            
            % change datatype for all 3D volumes to single precision float
            [v.dt] = deal([16 0]);
            
            % if wanted, change filename
            % [v.fname] = deal(niinewfile);
            
            % unset the private field, so no old file-descriptor is used!
            [v.private] = deal([]);
            
            % write volumes
            for vc = 1:numel(v)
                v(vc) = spm_write_vol(v(vc), y(:, :, :, vc));
            end
        end
        
    end
end
end

function [cellout] = stripemptycells(cellin)
if min(size(cellin))~=1 || max(size(size(cellin)))~=2
    error('Not a cell vector!')
end

keep = true(size(cellin));
for i = 1:length(cellin)
    if isempty(cellin{i})
        keep(i) = 0;
    end
end
cellout = cellin(keep);
end

function run_image_list = spm5_image_list(num_vols, file_list)
num_runs = length(num_vols);

if(length(file_list) == num_runs)
    % EITHER cell array with full image list OR single 4-D image names
    
    for i = 1:num_runs
        if size(file_list{i}, 1) == num_vols(i)
            % we have the full list already % tor edit april 2010
            for j = 1:num_vols(i)
                run_image_list{i}{j} = deblank(file_list{i}(j, :));
            end
        else
            % it's a single image name; expand it and create names
            % (could use expand_4d_images as well)
            printf_str = ['%' int2str(size(int2str(max(num_vols)), 2)) 'd'];
            for j = 1:num_vols(i)
                run_image_list{i}{j} = [deblank(file_list{i}) ',' num2str(j)];
            end
        end
    end
    
else
    % another format; not cell array of length(nruns)
    st = [0 cumsum(num_vols(1:end-1))];
    for i = 1:num_runs
        for j = 1:num_vols(i)
            run_image_list{i}{j} = [file_list{st(i) + j} ',1'];
        end
    end
end
end % function








