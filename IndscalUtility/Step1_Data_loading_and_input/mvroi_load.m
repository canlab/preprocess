function [DATA] = mvroi_load(varargin);
% function [DATA] = mvroi_load(varargin);
%
%% this routine tests whether there is a specifications file in your path, and if
%% there isn't, runs mvroi_specs so you can create one.
%% running mvroi with no inputs will always cal mvroi_specs
%% details of required specifications can be found by typing help mvroi_specs;
%
% e.g., DATA = mvroi_load('pain15')
% e.g., DATA = mvroi_load('pain15',pwd)

if length(varargin) == 0  % if there are no inputs to mvroi
    [DATA mypath]=mvroi_specs;   % make a new spec file
    prefix=DATA.SPEC.prefix;  
    
    if ~isfield(DATA.SPEC,'mypath'),DATA.SPEC.mypath = pwd;, end
    try
        % IT MIGHT BE BETTER NOT TO ADD TO THE PATH--COULD GET CONFUSING
        % W/MULTIPLE ANALYSES
        disp(['adding path ',DATA.SPEC.mypath]);
        eval(['addpath ',DATA.SPEC.mypath]);     %cd into the specified path                       
        eval(['cd ',mypath]);    %cd into the place specified    
        eval(['mkdir ',prefix]);  %make a dir with the prefix name
        eval(['cd ',prefix]);     %cd into there    
    catch
        warning('could not add path properly - is there a space character in your path?');
    end
elseif length(varargin) == 1    % if one input found
    prefix = varargin{1};
    if exist([prefix,'_specs.mat']);    % search the path for a spec file of that name
       % eval(['which ',prefix,'_specs.mat']);  %where did you find the spec file?
       % eval(['cd ',wh]);     %cd into there        
        disp('found specs....loading');
        load([prefix,'_specs']);
     
        if ~isfield(DATA.SPEC,'mypath'),DATA.SPEC.mypath = pwd;, end       
        try
            DATA.SPEC.prefix=prefix;
            disp(['adding path ',DATA.SPEC.mypath]);
            eval(['addpath ',DATA.SPEC.mypath]);     %cd into the specified path                
            eval(['cd ',DATA.SPEC.mypath]);     %cd into the specified path
            eval(['mkdir ',prefix]);  %make a dir with the prefix name
            eval(['cd ',prefix]);     %cd into there    
        catch
            warning('could not add path properly - is there a space character in your path?');
        end
    else 
        disp('couldnt find a spec file of that name...check path?');
        [DATA mypath]=mvroi_specs;   %run mvroi_specs
        prefix=DATA.SPEC.prefix;
      
        if ~isfield(DATA.SPEC,'mypath'),DATA.SPEC.mypath = pwd;, end      
        try
        disp(['adding path ',DATA.SPEC.mypath]);
        eval(['addpath ',DATA.SPEC.mypath]);     %cd into the specified path                
        eval(['cd ',mypath]);     % cd to the specified root path  
        eval(['mkdir ',prefix]);  % make a dir with the prefix name
        eval(['cd ',prefix]);     % cd into there
        catch
        warning('could not add path properly - is there a space character in your path?');
        end
    end
     
elseif length(varargin) == 2
    prefix = varargin{1};
    mypath = varargin{2};    % path hard coded
    eval(['cd ',mypath]);     % cd into there
    
    specname = fullfile(pwd,[prefix '_specs.mat']);
    if exist(specname); % if spec file found
        disp('found specs....loading');
        load(specname);
        try
            disp(['adding path ',DATA.SPEC.mypath]);        
            eval(['addpath ',DATA.SPEC.mypath]);     %cd into the specified path        
            eval(['cd ',DATA.SPEC.mypath]);     %cd into the specified path
            eval(['mkdir ',prefix]);  %make a dir with the prefix name
            eval(['cd ',prefix]);     %cd into there     
        catch
            warning('could not add path properly - is there a space character in your path?');
        end

    else 
        fprintf(1,'Couldn''t find a spec file called %s\n ...check path?\n',specname);
        [DATA mypath]=mvroi_specs;
        prefix=DATA.SPEC.prefix
        try
            disp(['adding path ',DATA.SPEC.mypath]);        
            eval(['addpath ',DATA.SPEC.mypath]);     %cd into the specified path        
            eval(['cd ',mypath]);     % cd to the specified root path  
            eval(['mkdir ',prefix]);  % make a dir with the prefix name
            eval(['cd ',prefix]);     % cd into there
        catch
            warning('could not add path properly - is there a space character in your path?');
        end
    end
end





% check for missing field names
% prompt for input if they're missing
% save SPEC file

 if ~isfield(DATA.SPEC,'HP'),DATA.SPEC.HP = input('Enter HP filter cutoff in s : ');, end
  if ~isfield(DATA.SPEC,'nruns'),DATA.SPEC.nruns = input('Enter number of runs : ');, end
   if ~isfield(DATA.SPEC,'TR'),DATA.SPEC.TR = input('Enter TR in s : ');, end
if ~isfield(DATA.SPEC,'scaling'),DATA.SPEC.scaling = [];, end

if exist('specname') == 1, disp('Saving SPEC file.'), eval(['save ' specname ' DATA']);,end



%%%% here we get the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data contained in the cell variable 'dat' in the format:
% SUBJECTS * {VOLUMES*REGIONS}
% the use of cell matrices allows variable numbers of volumes/subject
if ~isfield(DATA.SPEC,'datfilename'),DATA.SPEC.datfilename = prefix;, end

if exist([DATA.SPEC.datfilename,'_dat.mat']);  
disp('found dat....loading');
load([DATA.SPEC.datfilename,'_dat']);  
else    
    disp(['you need to have a file called ',DATA.SPEC.datfilename,'_dat.mat in your path']);
    disp('this will contain a cell variable called dat')
    error('which is a vector of SUBJECTS, each a cell containing {VOLUMES (images) * REGIONS} data')
end
DATA.dat=dat;
DATA.SPEC.numsub=length(DATA.dat);  %number of SUBJECTS
d=DATA.dat{1};
DATA.SPEC.numreg=size(d,2);     % number of REGIONS

%%%% here we get the onsets %%%%%%
% load data contained in the cell variable 'onsets' in the format:
% SUBJECTS * {BLOCKS*ONSETS}
% the use of cell matrices allows variable numbers of volumes/subject
if ~isfield(DATA.SPEC,'onsets');
if exist([DATA.SPEC.datfilename,'_onsets.mat']);  
disp('found onsets....loading');
load([DATA.SPEC.datfilename,'_onsets']);  
DATA.SPEC.onsets=onsets;
else    
    disp(['if you want the task to be removed, you need to have a file called ',DATA.SPEC.datfilename,'_onsets.mat in your path']);
    disp('this will contain a cell variable called onsets')
    disp('which is a vector of SUBJECTS * {BLOCKS * ONSETS}')
    warning('NO TASK WILL BE REMOVED!!!!');
    DATA.SPEC.onsets=[];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old code (now obsolete?): this converts the data from Chris' original format
%datorig = DATA.dat;
%dat = celldat2matrix(datorig);
%DATA.dat = [];
%for i = 1:size(datorig,1)   % for each subject
%    DATA.dat{i} = celldat2matrix(datorig(i,:));
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the names of each ROI
if exist([DATA.SPEC.datfilename,'_names.mat']);
disp('found names....loading');
load([DATA.SPEC.datfilename,'_names']);
else
    for r=1:DATA.SPEC.numreg;
        names{r}=(['r',num2str(r)]); % put in default names
    end
end
DATA.names=names;


% generate & plot DATA.SPEC.states, the vector assigning each VOLUME to a
% particular BLOCK
if isfield(DATA.SPEC,'states'),
    disp('Found existing DATA.SPEC.states coding for task states.')
else 
     disp('TRYING TO CREATE DATA.SPEC.states -- WARNING -- Check to see if this is right.')    

     for s=1:DATA.SPEC.numsub;
        clear states;
        d=DATA.dat{s};
        vols=1:size(d,1);
        for b=1:length(DATA.SPEC.spersess);
            states(DATA.SPEC.blockonsets(b):DATA.SPEC.blockonsets(b)+DATA.SPEC.spersess(b)-1)=DATA.SPEC.blocktypes(b);
            DATA.SPEC.states{s}=states';
        end
    end
end

% Save SPEC file again with states
if exist('specname') == 1, disp('Saving SPEC file.'), eval(['save ' specname ' DATA']);,end



doplot=1;
if doplot;
%%%% this plots a sample design, for subject 1 %%%%
figure;
s=DATA.SPEC.states{1};  % we only plot the first subject, as an example;
if ~isempty(DATA.SPEC.onsets);
    o=DATA.SPEC.onsets{1}; 
    for t=1:length(o);      %get event onsets for each trigger type
        oo(t,:)=o{t};
    end
    oo=oo(:);
    s(oo)=s(oo)+5;
end

plot(s);
set(gca,'Ylim',[min(s)-1 max(s)+1]);
Xlabel('scans');
title('example design (subject 1)');
try, saveas(gcf,'sample_design_matrix','tif');, close, catch,disp('Problem printing .tif of design matrix.'), end

% check states
try
    tmp = cat(2,DATA.SPEC.states{:}); 
    figure;imagesc(tmp), colormap gray, colorbar, xlabel('Subjects'),ylabel('Time'), title('Task States');,
    drawnow, pause(1)
    try, saveas(gcf,'task_states','tif');, close, catch, disp('Cannot save states image.'), end
catch, disp('Cannot display task states.  Are they different lengths for different subjects?'),
end

end


