function canlab_combine_structurals(basedir,struct_wildcard,structname)
%usage: canlab_combine_structurals(basedir,struct_wildcard)
%author: Scott Schafer
%date: 7/2/2010
%purpose:   This function takes multiple structural images taken at
%           different times, aligns them, and then creates an average that
%           is used for the coregistration of functional images.  The
%           output is designed for the format used in canlab_preproc.  
%           The mean structural image is saved in
%           basedir/[first run directory]
%
%           basedir: The base subject directory
%           struct_wildcard:  A wildcard string of where the structural
%               files to be combined are located.
%
%           If you want the mean file to move to basedir, add the
%           structname variable.  It should be the name of the mean file
%           that matlab creates. That name is 'mean' followed by the name
%           of the first structural image that you combine.
%           e.g., meanC04Wk0_struct.nii
%
%
%Example: if you're using a PC, you can only include one wildcard at the end, 
%               canlab_combine_structurals('/data/subject1/','Structs/Wk*/')
%             and matlab won't be able to find your images, so it'll open
%             the GUI
%         using a mac, you can include more than one wildcard
%               canlab_combine_structurals('/data/subject1/','Structs/Wk*/*.nii')
file_wildcard = fullfile(basedir, struct_wildcard);
structs = filenames(file_wildcard, 'char', 'absolute');
fprintf('Structural realignment started\n')

nimgs = size(structs,1);
matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = cell(nimgs,1);
for n = 1:nimgs
    matlabbatch{1}.spm.spatial.realign.estwrite.data{1}{n} = structs(n,:); %#ok
end
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

nrun = 1; % enter the number of runs here
% jobfile = {'Z:\Lane_ABRC\Controls\C03\structalign_job.m'};
% jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('serial', matlabbatch, '', inputs{:});

% spm_jobman('run',matlabbatch)

fprintf('Structural realignment finished \n')

if nargin > 2
    d = fileparts(structs(1,:));
    movefile(fullfile(d,structname),fullfile(basedir,structname))
end
