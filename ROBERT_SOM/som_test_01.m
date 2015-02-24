%
% A job to illustrate how the SOM data analysis works.
% 
% This assumes that the data is fMRI time series. 
%
% The job does the following:
%
%     1) gets a list of files to use in the analysis, all 
%        files determined by the directory specified by "theDIR", 
%        and with file root given by "theFILE". (See below).
%     2) using the list of files, P, a binary mask is created
%        with the same method as spm99. The file written is called
%        "mask.img" in the same directory that the data are found.
%     3) Data are read in, but only where the mask exists.
%     4) Data are detrended (linearly).
%     5) Data are high passed filtered, using parameters "highCutOff",
%        "fmriTR" and "nSessions". High pass filtering is the same as
%        found in spm99. 
%     6) SOM is calculated (see SOM_CalculateMap and other routines
%        on how this is done.)
%     7) Supercluster maps are made. These can be used for writing out
%        data or merging solutions etc.
%        

cp = cputime;
tic;

% Are the data broken into sessions.
nSessions  = 1;

% How long to iterate for.
nIter      = 50;

% Assuming a square NET, the length of the side.
nGrid      = 10;

% Sampling rate from the magnet.
fmriTR     = 2;

% High pass filter cutoff period.
highCutOff = 50;

theDIR = '/Users/rcwelsh/Experiments/TaylorLab/fcmriSimul/long_locked';
theFILE = 'fmri*.img';

% Grab the list of files to work upon.
P = spm_get('files',theDIR,theFILE);

cd(theDIR)

% Create the Mask on where to actually analyze data. Writes out mask.img
SOM_CreateMask(P);

[pn fn] = fileparts(P(1,:));
maskName = fullfile(pn,'mask.img');

% Read only the masked data. Returns an array of theData(nVoxels,nTime).
[theData maskInfo] = SOM_PrepData(P,'mask.img');

% Remove any linear trend - this also means it to zero.
theDataD = spm_detrend(theData',1)';

% Use the SPM routines for high pass filtering.
nScanPerSession = size(P,1)/nSessions;
for iSession = 1:nSessions
    K{iSession}.RT = fmriTR;
    K{iSession}.row = (iSession-1)*nScanPerSession+[1:nScanPerSession];
    K{iSession}.LChoice = 'none';
    K{iSession}.HChoice = 'specify';
    K{iSession}.HParam = highCutOff; % .02 Hz.
end

K = spm_filter('set',K);

fprintf('Filtering Data\n');

% High pass filter the data.
theDataF = spm_filter('apply',K,theDataD')';

% Now prepare to actually calculate the SOM.

nSOM = nGrid^2;

fprintf('Entering SOM Calculation\n');

% Calculate the SOM 

SOMResults        = SOM_CalculateMap(theDataF,nSOM,nIter);

% Store the header information - needed for writing out results as images.

SOMResults.header = maskInfo.header;
SOMResults.iMask  = maskInfo.iMask;

% Organize the data into super clusters.

[SOMResults.SuperCluster SOMResults.nCluster] = SOM_SuperCluster(SOMResults.SOM);

% Final amount of time.
toc
cputime - cp

%
% All done.
%