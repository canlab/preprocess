% convertPARtoNIFTI(runsDir, structFileIdentifier, PVEc, varargin)
% Converts ALL the .PAR and .par files found under a given directory to 3D nifti
% files.  Converts both functional and structural.  The program assumes
% that all files found are functional, UNLESS filename contains structFileIdentifier
%
% sample usage:
% convertPARtoNIFTI('/Volumes/RAID1/labdata/current/APP-fMRI', 'MPRAGE', 'pat1', 'WIP')
%   runsDir: searched for all .PAR and .par files under this, those files are
%   converted
%   struct identifier: unique pattern matching file name for struct file
%   varargin: series of parameters which are file patterns to skip (i.e.,
%   field maps)
%
% Skips any files in a directory which has _PAR_files as part of the name


function convertPARtoNIFTI(runsDir, structFileIdentifier, varargin)

    if ispc
        fprintf('This code only works on mac/unix currently, due to the call to the find command\n');
        exit
    end

     pfID = '_PAR_files';

    %append file sep to runsDir
    if runsDir(size(runsDir,2)) ~= filesep
        runsDir = [runsDir filesep];
    end

    %find all dirs which have .PAR and .par files
    command = ['find ' runsDir ' -iname \*.PAR'];
    [tmp,output] = system(command);
    
    %find all dirs which have .REC and .rec files
    command = ['find ' runsDir ' -name *.REC'];
    [tmp,recFiles] = system(command);   
    
    remain = output;
    dirlist = {};
    n = 1;
    while remain 
        [file, remain] = strtok(remain, char(10)); %split on newline
        %add file's dir to list of dirs
        temp=findstr('/',file);
        if temp
            current_path=fileparts(file);
            dirlist{n} = [current_path filesep];
            n = n+1;
        end
    end
    dirlist = char(dirlist);
    dirlist = unique(dirlist,'rows')

    excludePatterns = varargin;
    excludePatterns{end+1} = structFileIdentifier; %exclude struct as well
    
    for i=1:size(dirlist,1)
        currdir = deblank(dirlist(i,:))

        %skip par file dir
        if size(findstr(currdir, pfID)) > 0           
            continue
        end
        
        
        %convert structurals and functionals
        PAR_to_Nifti_BOLD_4pt2_dynsl_modified_v2_spm5(currdir, excludePatterns)
        PAR_to_Nifti_Structural_univ_4pt2_spm5(currdir,structFileIdentifier, 'n')
    end
  
    
    %move PAR and REC files aside
   parFiles = strrep(output, char(10), ' '); %list of .PAR files. replace newlines with spaces
   recFiles = strrep(recFiles, char(10), ' ');
   
  
   mkdir(runsDir, pfID);
   command = ['mv -v ' parFiles runsDir pfID]
   system(command);
   command = ['mv -v ' recFiles runsDir pfID]
   system(command);
 
     
    %tar and zip them .par files so don't get converted again later
    temp = findstr(runsDir, filesep);
    filename = [runsDir runsDir((temp(end-1)+1):end-1) pfID '.tar.gz']; %get a nice file name..
    command = ['tar -cvzf ' filename ' ' runsDir pfID filesep '*'];
    %system(command);
    
    fprintf('\n\nPlease confirm that %s contains all the PAR and REC files, then you may delete directory %s\n', filename, [runsDir 'PAR_files']);
    
end