% test for meta-analysis data.
% tor wager modified this from Robert Welsh's som_test_01.m



% Are the data broken into sessions.
nSessions  = 1;


theDIR = pwd;
theFILE = '*.img';


% Grab the list of files to work upon.
P = spm_get('files',theDIR,theFILE);

% Create the Mask on where to actually analyze data. Writes out mask.img
SOM_CreateMetaMask(P);

[pn fn] = fileparts(P(1,:));
%maskName = fullfile(pn,'mask.img');
maskName = spm_get(1,'*IMAGE','Select mask.');

% Read only the masked data. Returns an array of theData(nVoxels,nImages).
[theData maskInfo] = SOM_PrepData(P,maskName);

% normalize data



fprintf('Entering SOM Calculation\n');


cp = cputime;
tic;



% How long to iterate for.
nIter      = 50;

% Assuming a square NET, the length of the side.
% 

nGrid = round(sqrt(size(P,1) ./10));
disp(['nGrid: ' num2str(nGrid)]);



% Now prepare to actually calculate the SOM.

nSOM = nGrid^2;



% Calculate the SOM 

SOMResults        = SOM_CalculateMap(theData,nSOM,nIter);

% Store the header information - needed for writing out results as images.

SOMResults.header = maskInfo.header;
SOMResults.iMask  = maskInfo.iMask;

% Organize the data into super clusters.

[SOMResults.SuperCluster SOMResults.nCluster] = SOM_SuperCluster(SOMResults.SOM);

% Final amount of time.
toc
cputime - cp


mkdir SOM_Results
!mv mask* SOM_Results/
cd('SOM_Results');


% now make results masks
[wts,indices] = SOM_WriteIMGS(SOMResults);
wts = sort(wts',1,'descend');
tor_fig; plot(wts,'k','LineWidth',3); title('Map weights');
xlabel('Map')

SOM_ViewMap


%[S,H] = silhouette(theData, SOMResults.IDX,'Euclidean');

[cl,anyStudy,OUT] = meta_SOMclusters(SOMResults,theData,[],'group');

