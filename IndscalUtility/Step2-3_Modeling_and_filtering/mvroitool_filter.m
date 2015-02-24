function [M,cl] = mvroitool_filter(cl,varargin)
% [M,cl] = mvroitool_filter(cl,[images_per_session],[TR],[HPlength],[dummyimages])
%
% filtering and averaging of raw data within clusters
% for use before using matrix_multivariate
%
% returns M (3-D matrix input to matrix_multivariate)
% and cl structure with adjustedy field for each cluster.
%
% functional commands (inputs): Keyword followed by input variable
% case 'sessions', images_per_session
% case 'tr', TR 
% case 'hplength', HPlength 
% case 'dummy', dummyimgs = images at the start of each session todummy-code
% case 'trimstd', Windsorize at trimstd standard deviations
% 
% case 'inputfield', field name of cl to get raw data from
% case 'outputfield', field name of cl to save average data for each subject
%
% Two modes: single timeseries and cl
%
% Example for a single timeseries with dummy data:
% --------------------------------------------------
% images_per_session = [175   175   175   175   175   175];
% TR = 2;
% HPlength = 100;
% dummyimgs = 1:2;  
% yraw = randn(1050,1);
% y = mvroitool_filter(yraw,'sessions',images_per_session,'tr',TR,'hplength',HPlength,'dummy',dummyimgs)
%
% Example for clusters (cl = acc)
% % --------------------------------------------------

% ---------------------------------------------------------------------
% Set up inputs and defaults
% ---------------------------------------------------------------------

doplot = 0; verbose = 0; dummyimgs = []; trimstd = 3;
inputfield = 'raw_data'; outputfield = 'all_data';

for i = 1:length(varargin)
    if isstr(varargin{i})
        switch lower(varargin{i})
            % reserved keywords
            case 'betas'
            case 'design'

            % functional commands
            case 'sessions', images_per_session = varargin{i+1};
            case 'tr', TR = varargin{i+1};
            case 'hplength', HPlength = varargin{i+1};
            case 'dummy', dummyimgs = varargin{i+1};
            case 'trimstd', trimstd = varargin{i+1};
                
            case 'inputfield', dummyimgs = varargin{i+1};
            case 'outputfield', outputfield = varargin{i+1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isstruct(cl) && isfield(cl,inputfield)
    yraw = cl(1).(inputfield)(:,1,1);
else
    yraw = cl;  % data is input
end

% ---------------------------------------------------------------------
% Set up filter matrices for fast filtering
% if we enter only a timeseries as data, we're done
% ---------------------------------------------------------------------
[y,Xi,HP] = hpfilter(yraw,TR,HPlength,images_per_session,[],dummyimgs);
if ~isempty(trimstd) && trimstd > 0
    [y,OUT.ntrimmed,OUT.allw] = trimts(y,trimstd,[]);
end

if ~isstruct(cl)
    M = y;
    return
end

% ---------------------------------------------------------------------
% Filter and average data in clusters (in inputfield field)
% ---------------------------------------------------------------------
fprintf(1,'mvroitool_filter > ');

for i = 1:length(cl)
    
    fprintf(1,'Cl %03d, S %03d',i,0)
    
    % get sizes: standard format is time x voxels x subjects
    % (extract_raw_data.m)
    dat = cl(i).(inputfield);
    [t,v,s] = size(dat);
    
    % initialize average data output
    avg_data = zeros(t,s);
    
    for subj = 1:s        % subject
        
        fprintf(1,'\b\b\b%03d',subj);
        
        % initialize for this subject
        % remove voxels with empty data
        [rawdat,subjy,t,v] = get_subject_data(dat,subj);

        for voxel = 1:v   % voxel
            
            % filter this voxel
            y = hpfilter(rawdat(:,voxel),[],HP,images_per_session,Xi);
            [y,OUT.ntrimmed,OUT.allw] = trimts(y,trimstd,[]);
            
            subjy(:,voxel) = y;
        end
        
        submean = nanmean(subjy')';         % avg over voxels
        
        % add subject's mean to matrix of averages
        avg_data = add_subject(avg_data,submean,subj);
  
    end
    
    cl(i).(outputfield) = avg_data;
    fprintf(1,[repmat('\b',1,10) '%03d, S %03d'],i,0);
end

% get 3-D matrix form
for i = 1:length(cl), M(:,i,:) = cl(i).all_data;, end

return


% ---------------------------------------------------------------------
%
% Sub-functions
%
% ---------------------------------------------------------------------



function [rawdat,subjy,t,v] = get_subject_data(dat,subj)
rawdat = dat(:,:,subj);

% eliminate zero or empty voxels
whempty = sum(rawdat == 0 | rawdat == NaN);
whempty = find(whempty); % eliminate voxels with ANY 0 or NaN values
if any(whempty), fprintf(1,'\twarning: %3.0f empty voxels\n',length(whempty));,end
rawdat(:,whempty) = [];
if isempty(rawdat),
    fprintf(1,'\twarning: Subject has no valid voxels!!!\n');,
    subjy = [];
    return
end

subjy = zeros(size(rawdat));

[t,v] = size(subjy);

return




function avg_data = add_subject(avg_data,submean,subj);

[t,s] = size(avg_data);
[t2,v] = size(submean);

if v > 1, disp('Programming error: Problem in averaging process.'), keyboard;, end

d = t-t2;
if d > 0
    disp('Warning: subject''s timeseries too short.  Padding.');
    submean = padarray(submean,d,NaN,'post');
elseif d < 0
    disp('Warning: subject''s timeseries too long.  Truncating.');
    submean = submean(1:t,:);
end

avg_data(:,subj) = submean;

return
