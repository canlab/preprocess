function [imgs, rundirs, images_per_session, num_disdaqs] = canlab_preproc_list_files(basedir, run_wildcard, image_wildcard, disdaqs)
%[imgs, rundirs, images_per_session, num_disdaqs] = canlab_preproc_list_files(basedir, run_wildcard, image_wildcard, disdaqs)
%
% -------------------------------------------------------------------------
% List files and print
%
% Expand 4-D image filenames for SPM analysis if needed
% Check that files exist
% Move discarded acquisitions to subfolder, if requested
% -------------------------------------------------------------------------
%
% Used in canlab_preproc_2012.  See help canlab_preproc_2012 for info on
% inputs.


%% -------------------------------------------------------------------------
% List files and print
% -------------------------------------------------------------------------

str = 'CANlab preproc: Listing and checking files within each run.';
print_header1(str)

[imgs, rundirs, images_per_session, num_disdaqs] = deal([]);

numcards = size(run_wildcard, 1);

% establish each directory that has data

if numcards == 0  % No run wildcards; files are in raw directly
    error('You must place files in separate run subdirectories.');
    
elseif numcards == 1
    rundirs = dir(fullfile(basedir, 'Functional', 'Raw', run_wildcard));
elseif numcards == 2
    smalldirs = dir(fullfile(basedir, 'Functional', 'Raw', deblank(run_wildcard(1,:))));
    smallnum = numel(smalldirs);
    
    rundirs = [];
    for n=1:length(smalldirs)
        tempdirs = dir(fullfile(basedir, 'Functional', 'Raw', smalldirs(n).name, deblank(run_wildcard(2,:))));
        rundirs = [rundirs; dir(fullfile(basedir, 'Functional', 'Raw', smalldirs(n).name, deblank(run_wildcard(2,:))))]; %#ok
        
        smallernum = numel(tempdirs);
        if n==1, prevsmallernum = 0; end %doesn't matter if n=1 anyways
        
        for j = 1:smallernum
           ind = j+prevsmallernum*(n-1);
           rundirs(ind).name = fullfile(smalldirs(n).name,rundirs(ind).name);
        end
        prevsmallernum = smallernum;
    end
end


% eliminate selections that are not directories (i.e. files)
isnotdir = ~cat(1, rundirs.isdir);

if all(isnotdir)
    disp('Found no valid run directories.')
    searchin = fullfile(basedir, 'Functional', 'Raw', run_wildcard);
    fprintf('Searching in: %s\n', basedir);
    fprintf('Basedir: %s\nrun_wildcard: %s\n', basedir, run_wildcard);
    fprintf('Example of where files should be based on inputs:\n')
    fprintf('%s\n', fullfile(basedir, 'Functional', 'Raw', run_wildcard));
    return
end

rundirs(isnotdir) = [];


%ignore empty directories
ndirs = length(rundirs);
[temp_imgs, temp_wildcard] = deal(cell(ndirs, 1));
toremove = zeros(ndirs,1);
tokeep = zeros(ndirs,1);
for i = 1:ndirs
    temp_wildcard{i} = fullfile(basedir, 'Functional','Raw',rundirs(i).name, image_wildcard);
    temp_imgs{i} = filenames(temp_wildcard{i}, 'char', 'absolute');
    if size(temp_imgs{i},1) == 0
        toremove(i) = 1;
    else
        tokeep(i) = 1;
    end
end

toremove = find(toremove);
tokeep = find(tokeep);
fprintf('Found %3.0f run (session) directories:\n', numel(tokeep));
fprintf('Ignored %3.0f empty directories:\n', numel(toremove));

%calculate the number of separate runs
rundirs = rundirs(tokeep);
dirnames = char(rundirs(:).name);
nruns = length(rundirs);

%create image variables with placeholders
images_per_session = zeros(nruns, 1);
[imgs, img_wildcard] = deal(cell(nruns, 1));
str4d = cell(1, length(rundirs));

% -------------------------------------------------------------------------
% List files and remove disdaqs if called for
% -------------------------------------------------------------------------

for i = 1:nruns
    % get filenames
    %this is the image locations in wildcard format
    img_wildcard{i} = fullfile(basedir, 'Functional','Raw', rundirs(i).name, image_wildcard);
    
    % filenames is Wagerlab SCN toolbox function
    % create a char matrix with all image paths
    imgs{i} = filenames(img_wildcard{i}, 'char', 'absolute');
    
    if isempty(imgs{i})
        disp('Run %3.0f: Found no images!', i)
        fprintf('Looking in:\n%s\n', img_wildcard{i});
    end
    
    % If 4-D images, then print size of file and number of volumes (3-D
    % acquisitions)
    if size(imgs{i},1) == 1
        % assume 4-D
        str4d{i} = '4-D images';
        images_per_session(i) = scn_num_volumes(imgs{i});
    else
        str4d{i} = '3-D images';
        images_per_session(i) = size(imgs{i},1);
    end
    
    % Remove disdaqs from analyzed data and move to a separate folder
    if disdaqs
        %parse disdaq info
        
        if numel(disdaqs)==1
            num_disdaqs = disdaqs;
            num_vols = images_per_session(i) - disdaqs;
        elseif size(disdaqs) == size(images_per_session)
            num_disdaqs = disdaqs(i);
            num_vols = images_per_session(i) - disdaqs(i);
        else
            error('The size of disdaqs is not the same as the number of runs')
        end
        
        % remove acquisitions from 4-D or 3-D image sets
        if strcmp(str4d{i},'4-D images')  %remove 4D disdaqs
            d_imgs = remove_disdaq_vols(imgs(i),num_vols,num_disdaqs);
            [d f e] = fileparts(imgs{i});
            if ~exist(fullfile(d,'disdaqs'),'dir')
                mkdir(fullfile(d,'disdaqs'));
            end
            movefile(imgs{i},fullfile(d,'disdaqs',[f e]))
            imgs{i} = expand_4d_filenames(d_imgs);
            
        else  %remove 3D disdaqs
            d = fileparts(imgs{i}(1,:));
            if ~exist(fullfile(d,'disdaqs'),'dir')
                mkdir(fullfile(d,'disdaqs'));
            end
            
            for n = 1:num_disdaqs %assume they're in order, and the top ones are earliest.
                [d f e] = fileparts(imgs{i}(n,:));
                movefile(imgs{i}(n,:),fullfile(d,'disdaqs',[f e]));
            end
            
            %update imgs to exclude disdaqs
            imgs{i} = imgs{i}((num_disdaqs+1):end,:);
            
        end
    else
        num_disdaqs = 0;
        if strcmp(str4d{i},'4-D images')
            imgs{i} = expand_4d_filenames(imgs{i});
        end
    end
    images_per_session(i) = size(imgs{i}, 1);
    
    fprintf('Run %3.0f:\t%s\t%3.0f images\t%3.0f discarded\t%s\n', i, str4d{i}, images_per_session(i), num_disdaqs, img_wildcard{i});
    
end % nruns

end % function



function print_header1(str)

s = '======================================================================================';
len = length(s);

disp('======================================================================================');
disp('=                                                                                    =');
fprintf('= %s%s=\n', str, repmat(' ', 1, max([1 length(s) - length(str) - 3])));
disp('=                                                                                    =');
disp('======================================================================================');


end
