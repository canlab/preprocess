function [allerrors, ISDONE] = canlab_preproc_check_PREPROC(PREPROC, varargin)
% allerrors = canlab_preproc_check_PREPROC(PREPROC, [optional inputs])
%
% Checks for valid fields and files in PREPROC object, used in
% canlab_preproc.m
%
% You can enter a field, a description for the field, and a flag (1/0) for
% "error_on_missing" (for files)
% Example: check for ra_images before normalization:
% allerrors = canlab_preproc_check_PREPROC(PREPROC, 'ra_func_files', 'Realigned images to warp', 1)
%
% Or, by default, it checks these fields, and returns errors for fields
% that should not be empty missing at initial setup (marked with *).
%
% optional input: 'erroroff', it only checks these fileds and doesn't
% return errors. (added by Wani 09/01/12)
% Example: [errors, ISDONE] = canlab_preproc_check_PREPROC(PREPROC, 'erroroff');
%
% PREPROC.meanfilename
% PREPROC.meanfilename_realigned
% * PREPROC.func_files
% * PREPROC.images_per_session
% PREPROC.implicit_mask_file
% * PREPROC.anat_files
% PREPROC.ra_func_files
% PREPROC.wra_func_files
% PREPROC.swra_func_files

allerrors = {};
ISDONE = [];

fields = {'anat_files' 'func_files' 'implicit_mask_file' 'meanfilename' 'meanfilename_realigned' 'ra_func_files' 'wra_func_files' 'swra_func_files'};
descrip = {'Anatomical' 'Functional before preprocessing' 'Implicit mask before preprocessing' 'Mean before realign' 'Mean after realign' 'Realigned functionals' 'Warped/realigned functionals' 'Smoothed/warped/realigned functionals'};
    
if nargin < 2
    error_on_missing = [1 1 0 0 0 0 0 0 0];
    
elseif nargin == 2 && strcmp(varargin{1}, 'erroroff')
    error_on_missing = [0 0 0 0 0 0 0 0 0];
    
else
    fields = varargin(1);
    descrip = varargin(2);
    error_on_missing = varargin{3};
end


for i = 1:length(fields)
    
    if ~isfield(PREPROC, fields{i})
        PREPROC.(fields{i}) = [];
        ISDONE.(fields{i}) = 0;
    end
    
    myvar = PREPROC.(fields{i});
    
    [currenterrors, ISDONE.(fields{i})] = preproc_check_field(myvar, fields{i}, descrip{i}, error_on_missing(i));
    
    allerrors = [allerrors; currenterrors];
    
end % fields

if ~isempty(allerrors)
    print_header1('Errors found in PREPROC setup. Quitting.');
    disp(char(allerrors{:}))
    error('Fix PREPROC before running.');
end

end % function










function [errors, isdone] = preproc_check_field(myvar, fieldname, descrip, error_on_missing)

errors = {};
dispstr = '';
isdone = 0;

switch fieldname
    
    case {'anat_files' 'meanfilename' 'meanfilename_realigned' 'implicit_mask_file'}
        
        if iscell(myvar)
            errors{end + 1, 1} = sprintf('PREPROC.%s should not be a cell array.\n', fieldname);
            myvar = cellstr(myvar); % for later checking
            ismissing = 0;

        elseif exist(myvar, 'file')
            dispstr = sprintf('PREPROC.%s\t%s\n\tFile exists.\t\n', fieldname, descrip);
            ismissing = 0;
            isdone = 1;
            
        elseif isempty(myvar)
            dispstr = sprintf('PREPROC.%s\t%s\n\tDoes not exist yet.\t\n', fieldname, descrip);
            ismissing = 1;
             
        else % Missing
             dispstr = sprintf('PREPROC.%s\t%s\n\tIs named in PREPROC but does not exist.\t\n', fieldname, descrip);
            ismissing = 1;
        end
        
        fprintf(dispstr)
        
        if ismissing && error_on_missing
            errors{end + 1, 1} = sprintf('Missing required file in PREPROC.%s\n', fieldname);
        end
        
        
    case {'func_files' 'a_func_files' 'ra_func_files' 'wra_func_files' 'swra_func_files'}
        
        if isempty(myvar)
            dispstr = sprintf('PREPROC.%s\t%s\t\n\tDoes not exist yet.\n', fieldname, descrip);
            ismissing = 1;
            fprintf(dispstr)
            
            if ismissing && error_on_missing
                errors{end + 1, 1} = sprintf('Missing required file (empty field) in PREPROC.%s\n', fieldname);
            end
            
        elseif ~iscell(myvar)
            errors{end + 1, 1} = sprintf('PREPROC.%s must be a cell array with expanded image names.\n', fieldname);
            ismissing = 0;
            
        else
            dispstr = sprintf('PREPROC.%s\t%s\n', fieldname, descrip);
            fprintf(dispstr)
            
            for i = 1:length(myvar)
                myfilename = getfrom4dnames(myvar{i});
                ismissing(i, 1) = ~exist(myfilename, 'file');
                
                if ~ismissing(i, 1)
                    dispstr = sprintf('\tFile exists.\t%s\n', myfilename);
                else
                    dispstr = sprintf('\tDoes not exist yet.\t%s\n', myfilename);
                end
                
                fprintf(dispstr)
                
                if ismissing(i) && error_on_missing
                    errors{end + 1, 1} = sprintf('Missing required file in PREPROC.%s\n', fieldname);
                end
                
            end % each cell
            
            if ~any(ismissing), isdone = 1; end
            
        end
        
        
        
        
    otherwise
        error('canlab_preproc_check_PREPROC: Unknown field %s\n', fieldname)
end

end% function



function  myfilename = getfrom4dnames(myimages)

n = size(myimages, 1);
allframes = cell(n, 1);

for i = 1:n % a string matrix of file names, potentially with ,x for SPM frame no.
    
    myframe = remove_image_frame_comma(myimages(i, :));

    allframes{i, 1} = myframe;
    
end

allframes = char(allframes{:});
myfilename = unique(allframes, 'rows'); % unique files

end


function myfile = remove_image_frame_comma(myfile)

wh = find(myfile == ',');
    
    if any(wh)
        wh = wh(1);
        myfile(wh:end) = [];
    end


end




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


