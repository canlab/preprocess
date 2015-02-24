function [CLU,clusters,subclusters] = snpm2struct(cwd,varargin)
% [CLU,clusters,subclusters] = snpm2struct(cwd, [opt] reslice, [opt] alt. P)
%
% Given the directory of an SnPM analysis,
% Creates a structure with relevant info
% Compatible with spm_add_blobs and cluster functions.
%
% to display:
% spm_image	(and choose anatomical) or spm_image('init','myname.img')
% spm_orthviews('AddBlobs',1,CLU.XYZ,CLU.Z,CLU.V.mat)
% spm_orthviews('AddColouredBlobs',1,CLU.XYZ,CLU.Z,CLU.V.mat,[0 0 1])
%
% optional: reslice flag (enter 0, 1, or 2)
% if selected, this option will allow you to 
% reslice the SnPM filtered t into another space
% and choose new images to extract data from.
% 
% entering 1 will choose new images only,
% and entering 2 will perform the reslice and get new images.
% 
% Tor Wager, 2/12/02

if length(varargin) > 0, resl = varargin{1};,else resl = 0;,end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% go to the input directory and check for results files
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

eval(['cd(''' cwd ''')'])

% this is the file to load with the filtered results
P = 'SnPMt_filtered.img';

% reslice if necessary
if resl == 2, reslice_imgs([],P,0);, P = ['r' P];,end

if ~(exist(P) == 2), error(['Cannot find ' P ' in directory.']), end
if ~(exist('SnPM.mat') == 2), error('Cannot find SnPM.mat file.'), end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% get coordinates and height values from the mask
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

[CLU.XYZ,CLU.XYZmm,CLU.Z,CLU.V] = img2voxel(P);
CLU.XYZ = CLU.XYZ(1:3,:);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% change pseudo-T values to Z scores
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load SnPM                       % to get df
[CLU.Z] = spm_t2z(CLU.Z,df);


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% fill in other fields
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% adjust crit_t from 0 to 1 for compatibility with fixed_TSU (Talairach)
% needs to be a non-zero value to display on Talairach
% ------------------------------------------------------
CLU.crit_t = ST_Ut;

CLU.cl_size = 0;
CLU.df = df;

% for compatibility with SPM struct and cluster analysis
% ------------------------------------------------------
CLU.voxSize = diag(CLU.V.mat)'; 
CLU.voxSize = CLU.voxSize(1:3);
CLU.VOX = CLU.voxSize;		% compatible with VOL structure
CLU.M = CLU.V(1).mat;

CLU.u = CLU.crit_t;
CLU.k = CLU.cl_size;

[a,b,c] = fileparts(CLU.V.fname);
CLU.title = ['snpm_' b];

if resl, P = spm_get(Inf,'*img','Select imgs for data extraction.');
else, load SnPMcfg                % to get P
end

try,
    vv = spm_vol(P(1,:));
catch,
    warning('Image files cannot be found in ind. subject directories: Path has changed?')
    % tries to re-find P - only works for my specific directory 
    % takes input from EXPT file!
    P = varargin{2};
    disp('Taking P from EXPT input.')
    %for i = 1:size(P,1)
    %    try,
    %        vv = spm_vol(P(i,:));
    %        newP{i,1}=P(i,:);
    %    catch
   %         [d,f,e] = fileparts(P(i,:));
   %         [dummy,d] = fileparts(d);
   %         newP{i,1}=fullfile('..',d,[f e]);
   %     end
   %end
    %P = cell2mat(newP);
    %disp(['Found: ' P(1,:) ' etc.'])
    P
end

[clusters] = tor_extract_rois(P,CLU,CLU,1);
for i = 1:length(clusters), clusters(i).P = P;, end

[subclusters] = tor_get_spheres3(clusters);
return