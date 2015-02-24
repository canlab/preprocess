function canlab_spm_delete_betas(modeldirs, varargin)
% canlab_spm_delete_betas(analysis_directories [options])
%
% DESCRIPTION
%   Delete beta images (ANALYZE and or NIfTI) in SPM subject level
%   analyses based on their name (given by SPM.xX.name).
%
%   Names are picked out with a regular expression (help regexp). The 
%   default expression is meant to pick out nuisance regressors (regressors 
%   loaded by the Multiple Regressors method). Change the expression to 
%   delete other sets of betas.
%
%   By default, will merely print the rm commands constructed.
%   You should at least start it this way first to make sure 
%   it's doing what you want, only then use the 'delete' option.
%
%   When done running, will print out number of files targeted and unique
%   names of files targeted. "Targeted" here doesn't mean the file did exist, 
%   but that it was deleted if it did.
%
% ARGUMENTS
%   analysis_directories
%       cell array of filenames (absolute paths preferred) of
%       SPM lower level analyses.
%
% OPTIONS
%   'e', regexp
%       pick out regressors matching regexp
%       DEFAULT: ' R[0-9][0-9]*$'
%   'v'
%       verbose mode: print individual rm commands
%       do still print summary upon completion (see above)
%   'delete'
%       actually delete image files
%       DEFAULT: only print out what would be done
%
% EXAMPLES
%   canlab_spm_delete_betas({'analysis/first_level/ilcp19'})
%   canlab_spm_delete_betas({'analysis/first_level/ilcp19'}, 'delete')
% or
%   sublevs = filenames('analysis/first_level/ilcp[0-9]*', 'absolute');
%   canlab_spm_delete_betas(sublevs, 'e', 'NLR_[^h]*')
%   canlab_spm_delete_betas(sublevs, 'e', 'NLR_[^h]*', 'delete')
%


if ~isunix
    fprintf('Sorry, this will only work on a linuxy system.\n')
    return   
end

%% setup
STARTINGDIR = pwd;

% defaults
DELETE = false;
VERBOSE = false;
expr = ' R[0-9][0-9]*$';

% make modeldirs are absolute
if ~exist('modeldirs','var')
    error('No SPM subject level directories specified.')
end

% parse other arguments
i=1;
while i<=numel(varargin)
   if ischar(varargin{i})
       switch varargin{i}
           case 'regexp'
               i=i+1;
               expr = varargin{i};
           case 'v'
               VERBOSE = true;
           case 'delete'
               DELETE = true;
           otherwise
               error(['UNRECOGNIZED OPTION: ' varargin{i}])
       end
   else
       disp(varargin{i})
       error('Above option UNRECOGNIZED')
   end
   i=i+1;
end


%% main
count=0;
for i = 1:numel(modeldirs)
    clear SPM regs
    try
        cd(modeldirs{i})
        load SPM.mat
        fprintf('MOVING to: %s\n',modeldirs{i})        
    catch exc
        disp(getReport(exc,'extended'))
        fprintf('SKIPPING PROBLEMATIC DIRECTORY: %s\n',modeldirs{i})
        cd(STARTINGDIR)
        continue
    end
        
    regs = find(~cellfun('isempty',regexp(SPM.xX.name,expr)));
    targfiles{i} = cell(numel(regs),1);    
    for r=1:numel(regs)
        targfiles{i}{r} = SPM.xX.name{regs(r)};
        cmd = sprintf('!rm beta_%04d.{img,hdr,nii} 2>/dev/null  # deleting %s',regs(r),SPM.xX.name{regs(r)});
        if VERBOSE, disp(cmd), end
	if DELETE, eval(cmd), end
        count=count+1;
    end
    cd(STARTINGDIR)
end


%% clean up
fprintf('\n\nNumber of betas targeted: %d',count)
allfiles = {};
for i=1:numel(targfiles)
    allfiles = [allfiles;targfiles{i}];
end
fprintf('\nNames of betas targeted:\n')
disp(unique(allfiles))

cd(STARTINGDIR)
