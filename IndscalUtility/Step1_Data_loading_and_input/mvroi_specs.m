function [SPEC]=mvroi_specs(SPEC,varargin)
% function [SPEC]=mvroi_specs(SPEC,[DATA])
% in this script we load defaults for mvroi analysis.
% each option is commented below
%
% MVROI analysis tool is for explore patterns of functional connectivity
% between pre-defined regions of interest (ROIs) for block design fMRI
% experiments.  The data to be input (matlab variable dat) consists of a cell vector,
% in which each cell is a matrix of VOLUMES * REGIONS for a single SUBJECT.
% this should be saved in your path as a .mat file named 'prefix_dat' where
% prefix is the name of your experiment.
%
%
% some terminology for BLOCKS SESSIONS and ONSETS:
%
% SESSION: We refer to the SESSION as the number of runs (scanner on/images acquired/scanner off) in 
% which the data were acquired.  
% 
% BLOCK: A BLOCK is defined by a separate vector for each subject (in SPEC.states) 
% which  assigns each volume to one or other condition (this will be generated automatically 
% from the inputs entered).  A typical SESSION thus consists of 2
% or more BLOCKS which may be cycled many times across the session or
% follow each other consecutively.
%
% ONSET: Each block will consist of several ONSETS, 
% which are specified timings of stimulus events which occurred during the 
% BLOCK.  BLOCKs typically differ with regard to
% a variable of interest and results will reflect how connectivity between
% REGIONS differs as a function of block.
% Event-related designs, in which different trial types follow each other
% in rapid random order, are not suitable for analysis with this toolbox.
%
%
% some terminology for REGIONS CLUSTERS and DIMENSIONS 
%
% if 2nd argument, this is DATA, and we get a few more things we need (at
% end)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('welcome to mvroi analysis toolbox.  type help mvroi_specs for more info');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('...');
disp('...');
disp('...');



disp('...');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('BLOCKS,ONSETS,COMPARISONS');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('...');

if ~isfield(SPEC,'nruns'), SPEC.nruns=input('how many experimental SESSIONS were there?...(eg 2)...');, end
% A SESSION is defined as a single run (scanner on; images acquired;
% scanner off).  An experiment may have multiple sessions.  It is necessary
% to model each session separately in order to remove scanner effects such
% as drift during preprocessing.
% we assume that sessions are of equal length (in VOLUMES).

if ~isfield(SPEC,'spersess'),SPEC.spersess = input('how many VOLUMES in each SESSION (or SCANNER RUN)?...(eg [100 100])....'); ,end
% Number of images in each SESSION, [n1 n2 ... nx]

% States
% --------------------------------
if ~isfield(SPEC,'states'),
    fprintf(1,'No SPEC.states field entered.  This is used to compare cov btwn blocks, states, conditions of interest.\n')
    fprintf(1,'states is a cell array, on cell per subject, 1 column vector per subject of time x 1 with integer codes for state/block types.\n')
    statefile = input('Do you want to load one?  Type the name of the file (in current dir) with states, no .mat extension: ','s');
    while ~exist(statefile) ==2 & ~isempty(statefile)
        statefile = input('Can''t find that file.  Try again, or return to continue: ','s');
    end
    if ~isempty(statefile), load(statefile), try, SPEC.states = states;, catch, error('You need a var called states in this file.'), end, end
end
 
if isfield(SPEC,'states'), 
    SPEC.numstates = max(cat(1,SPEC.states{:})); , 

    % check states
    try
        tmp = cat(2,SPEC.states{:}); 
        figure;imagesc(tmp), colormap gray, colorbar, xlabel('Subjects'),ylabel('Time'), title('Task States');,
        drawnow, pause(1)
        try, saveas(gcf,'task_states','tif');, close, catch, disp('Cannot save states image.'), end
    catch, disp('Cannot display task states.  Are they different lengths for different subjects?'),
    end
 
    if ~isfield(SPEC,'tasknames'),
        for c=1:SPEC.numstates
            SPEC.tasknames{c} = input(['what should I call condition ',num2str(c),'?...e.g. control...'],'s');
        end
    end
else
    fprintf(1,'One state only entered.\n')
    SPEC.numstates = 1;
end


% Global scaling and cross-correlation options
% --------------------------------
if ~isfield(SPEC,'doglobal'), SPEC.doglobal = input('Remove global mean? Type regression or none: ','s');,end
  
if ~isfield(SPEC,'betaflag'), SPEC.betaflag = input('Use betas instead of correlation z'' scores (rec. 0), 1/0: ');,end
    
if ~isfield(SPEC,'robustflag'), SPEC.robustflag = input('Use robust IRLS fitting (rec. 1), 1/0: ');,end
    

disp('saving mvroi_specs'); save mvroi_specs SPEC

% Onsets
% --------------------------------
if ~isfield(SPEC,'onsets'),
    fprintf(1,'\nNo SPEC.onsets field entered.  This is used to define where events of interest happen.\n')
    fprintf(1,'onsets is a cell vector, 1 cell per person, each cell contains a separate cell for each session, containing the onsets for each event type.\n')
    onsetfile = input('Do you want to load one?  Type the name of the file (in current dir) with onsets, no .mat extension: ','s');
    while ~exist(onsetfile) ==2 & ~isempty(onsetfile)
        onsetfile = input('Can''t find that file.  Try again, or return to continue: ','s');
    end
    if ~isempty(onsetfile), load(onsetfile), try, SPEC.onsets = onsets;, catch, error('You need a var called onsets in this file.'), end, end
end

if ~isfield(SPEC,'comps'),SPEC.comps=input('which comparison (contrast) weights?...(eg [-1 1])...');  ,end
while length(SPEC.comps) ~= SPEC.numstates
    SPEC.comps = input('SPEC.comps is incorrect. You need one contrast weight for each state.  Try again: ');
end

% Behavior
% --------------------------------
if ~isfield(SPEC,'beh'), 
    fprintf(1,'\nNo SPEC.beh field entered.  This is used behavioral predictors that may predict network covariance.\n')
    fprintf(1,'It is a column vector with one score per subject, called beh in a mat file, in the same order as in the data (dat).\n')
    behfile = input('Do you want to load one?  Type the name of the file (in current dir) with states, no .mat extension: ','s');
    while ~exist(behfile) ==2 & ~isempty(behfile)
        behfile = input('Can''t find that file.  Try again, or return to continue: ','s');
    end
    if ~isempty(behfile), load(behfile), try, SPEC.beh = beh;, catch, error('You need a var called beh in this file.'), end, end
end


% this defines contrasts to test connectivity between conditions.
% you can have as many contrasts as you like, but there should be an
% equal number of weights as conditions and they should sum to 1.
if ~isfield(SPEC,'comptitle'),SPEC.comptitle = input(['what should I call this comparison?...eg control>test...'],'s');,end
if iscell(SPEC.comptitle), SPEC.comptitle = SPEC.comptitle{1};,end

if ~isfield(SPEC,'firpoints'),SPEC.firpoints = input('vector of beta lengths with which to model each ONSET..(eg [16 16 16....16])....');,end
% vector or integer of number of betas to estimate for each onset across the entire session 
%(in TRs); same as EXPT.DX.numframes
% note that for most (but not all) experiments this should be ~32s/TR
% also note that this will be treated identically across all subjects


% load the names of each ROI
% --------------------------------
if ~isfield(SPEC,'names'),
    fprintf(1,'No SPEC.names field entered.  This contains short names for each region, cell array, 1 cell per region.\n')
    namesfile = input('Type the name of the file (in current dir), no .mat extension, or return to skip: ','s');
    while ~exist(namesfile) ==2 & ~isempty(namesfile)
        namesfile = input('Can''t find that file.  Try again, or return to continue: ','s');
    end
    if ~isempty(namesfile), load(namesfile), try, SPEC.names = names;, catch, error('You need a var called names in this file.'), end, end
end
    
if ~isfield(SPEC,'names'),
    for r=1:100;
        names{r}=(['r',num2str(r)]); % put in default names
    end
    SPEC.names=names;
end

% check names
for i = 1:length(SPEC.names), if isempty(SPEC.names{i}), SPEC.names{i} = ['r' num2str(i)];,end, end
    

disp('saving mvroi_specs'); save mvroi_specs SPEC



disp('...');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('DIMENSIONALITY, CORRELATION AND CLUSTERING');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('...');

if ~isfield(SPEC,'ndims'),SPEC.ndims={input('how many DIMENSIONS to the data?...(eg 3 or choose)...','s')};, end       
% number of dimensions for indscal solution
% enter an integer or 'choose' (string) to use PCA to choose ndims
% note that solutions of >3 dimensions can only be plotted as 3D.

if ~isfield(SPEC,'clustsolutions'),SPEC.clustsolutions=input('what range of CLUSTERS would you like?...(eg 2:10 or choose)...','s');, end
%how many clusters is reasonable?
% enter range of integers (3:12) or 'choose' (string) to automate: 
% for r  regions chooses 2:(r ./ 3) as max.

if ~isfield(SPEC,'numperm'),SPEC.numperm=input('how many permutations for testing?...(eg 5000)...');  ,end          
% number of permutations in cluster testing for boostrapping
% 1000 will invariably be the best option

if ~isfield(SPEC,'shiftby'),SPEC.shiftby=input('how many TRs to shift when cross-correlating?...(eg 2)...'); ,end    
% how many scans in each direction for sliding correlation?
% 0 gives no latitude for shifting; enter integer in TRs.
                         
if ~isfield(SPEC,'corrscaling'),SPEC.corrscaling={input('how should I scale the data?...eg rect or scalar or just press return to skip this...','s')};, end
% this can be set to 'rect' (rectify negative correlations) or 'scalar'
% (scale all correlations 0-1)

if ~isfield(SPEC,'remove'),SPEC.remove={input('what should I do about ELEMENTS which fit the CLUSTER poorly?...eg keep/iter/thresh....','s')};,end
% this is an input to testclust:
% 'keep' - testcluster.m, assumes that all ELEMENTS are included in the
% solution, and calculates cluster quality accordingly
% 'iter' - removes ELEMENTS which are below 5% of the permuted
% distribution, and recomputes quality.  NOTE that p-values for this
% solution will not be meaningful, as you are 'pruning' the solution to be
% significant by removing the poor fitting ELEMENTS
% 'thresh' removes ELEMENTS which fail to achieve 95% confidence for the
% permuted solution

disp('saving mvroi_specs'); save mvroi_specs SPEC


disp('...');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('FILTERING ETC');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('...');

if ~isfield(SPEC,'trim'),SPEC.trim=input('z-score for trimming and spike removal...(eg 3)...');   ,end     
% number of sd's to trim during preprocessing (for spike removal)

if ~isfield(SPEC,'TR'),SPEC.TR=input('volume sampling?  [TR]....(eg 2)...');, end
% the TR

if ~isfield(SPEC,'HP'),SPEC.HP=input('high pass filter in seconds....(eg 120)...');, end
% what high pass filter?

disp('saving mvroi_specs'); save mvroi_specs SPEC

% DATA-specific arguments
if length(varargin) > 0,
    DATA = varargin{1};
    
    SPEC.numsub = size(DATA.DATA(1).dat,2);
    SPEC.numreg = size(DATA.DATA(1).dat,1);
    try, SPEC.names = SPEC.names(1:SPEC.numreg);, catch, error('Names vector seems to be too short!'),end

    if ~isfield(SPEC,'states')
        fprintf(1,'Creating States: Assuming single task condition.\n');,
        SPEC.states{1} = ones(size(DATA.DATA.dat{1},1),1);,
    end
    
    if length(SPEC.states) < length(DATA.DATA(1).dat),
        fprintf(1,'Fewer states entries than subjects: replicating 1st subject.\n');,
        for i = 1:length(DATA.DATA(1).dat), SPEC.states{i} = SPEC.states{1};,end
    end
end



disp('...');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('SAVE SPECS');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('...');

%saver=input('save these specs in DATA????  warning: may overwrite...(eg y)...','s');
%if saver=='y';
    disp('Saving mvroi_specs.mat')
    save mvroi_specs SPEC
%else
%    disp('nothing saved')
%end


