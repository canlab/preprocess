function EXPT = wb_hewma(EXPT,varargin)
% EXPT = wb_hewma(EXPT,[which cons],[dools],[mask])
%
% Tor Wager, 10/2/05
%
% This function does a robust fit on contrasts specified in EXPT.SNPM.connames,
% using images in EXPT.SNPM.P
% The model should be stored (without intercept) in EXPT.cov
% Empty EXPT.cov will result in a one-sample t-test
%
% For each contrast, it saves con, t, p, filtered t imgs, and clusters.mat 
% p-value images are 2-tailed!
%
% All predictors should be in a matrix of column vectors, stored in EXPT.cov.
% They are centered and scaled in this function.
% Intercept is added as 1st predictor
%
% see get_expt_info.m for how to make the EXPT structure
% fields should be:
% EXPT.SNPM.P           containing image names of individual subjects
% EXPT.SNPM.connames    str mtx of contrast names for each P string matrix
% EXPT.SNPM.connums     contrast numbers for each P
% EXPT.cov              empty for 1-sample ttest, or containing covariates
% [mask] is the name of the mask file to use
%
% [Which cons] optional vector of which arrays in EXPT.SNPM.P to use
% e.g., [1 2 4] runs the first, second, and 4th sets of image names
%
% Example:
% EXPT = robfit(EXPT,[6:2:29]);   
%
% DO OLS, cons 1-4, use gray matter mask
% mask = which('scalped_avg152T1_graymatter_smoothed.img');
% [dummy,mask] = reslice_imgs(EXPT.SNPM.P{1}(1,:),mask);   % reslice mask
% to image space
% EXPT = robfit(EXPT,1:4,1,mask);

global xdim
global ydim
global emptyimage
global V
global domask
global lam
global f

fprintf(1,'-------------------------------\n')
fprintf(1,'Hierarchical Ewma - whole brain\n')
fprintf(1,'-------------------------------\n')

% -----------------------------------------------------
% Get sample image to get dimensions
% -----------------------------------------------------
if exist(EXPT.im_files{1}(1,:)) == 2
    V = spm_vol(deblank(EXPT.im_files{1}(1,:)));
    xdim = V.dim(1);
    ydim = V.dim(2);
    zdim = V.dim(3);
    emptyimage = zeros([xdim ydim zdim]);       % for writing slices of stat values
    disp(['Found sample image: ' V.fname]);
else
    error('Enter an image file in EXPT.im_files{1}(1,:) to get dimensions of data images.');
end
   
%load lambda
try
    load([EXPT.subjects{1} filesep 'lam.mat']);
    disp(['Lambda found: lam = ' num2str(lam)]);
catch
    disp(['Cannot find lambda parameter saved in lam.mat in ' EXPT.subjects{1}])
end

% -----------------------------------------------------
% Get slice .mat file names
% -----------------------------------------------------
if ~isfield(EXPT,'FILES'), EXPT.FILES = [];,end
if ~isfield(EXPT.FILES,'ewma_z'), 
    disp('Cannot find EXPT.FILES.ewma_z.  Looking for z_slice*.mat');
    EXPT = getfunctnames2(EXPT,'z_slice*.mat','FILES.ewma_z');  
else
    disp('Found EXPT.FILES.ewma_z.  Assuming z_slice*.mat images are there and using those.');
end

if ~isfield(EXPT.FILES,'ewma_var'), 
    disp('Cannot find EXPT.FILES.ewma_z.  Looking for var_*.mat');
    EXPT = getfunctnames2(EXPT,'var_*.mat','FILES.ewma_var');
else
    disp('Found EXPT.FILES.ewma_var.  Assuming var_*.mat images are there and using those.');
end

if zdim ~= size(EXPT.FILES.ewma_z{1},1)
    disp('Error: Number of slices in EXPT.FILES.ewma_z{1} does not match that of sample image file.');
end


% -----------------------------------------------------
% center and scale predictors
% -----------------------------------------------------
if ~isfield(EXPT,'cov'), EXPT.cov = [];, end
covt = EXPT.cov;

if ~isempty(covt)
    disp('Found behavioral covariate.')
    if length(covt) > size(covt,1), covt = covt';, end
    
    covt = scale(covt);
    
    %covt = covt - repmat(mean(covt),size(covt,1),1);
    %for i = 1:size(covt,2)
    %    covt(:,i) = covt(:,i) ./ var(covt(:,i));
    %end
else
    disp('<no behavioral covariates>')
end

% intercept: automatically added if doing robustfit.m
% covt = [ones(1,size(covt,1)) covt];
%covt(:,end+1) = 1;

EXPT.cov = covt;


% ----------------------------------------------------
% * process optional input args
% ----------------------------------------------------

if length(varargin)>0,wh=varargin{1};,else,wh=1:length(EXPT.FILES.ewma_z);,end
if length(varargin)>1,dools = varargin{2};, else, dools = 1;, end
if length(varargin)>2,domask = varargin{3};, else, domask = 0;, end

if domask,disp(['Using mask: ' domask]), else, disp('<no mask image specified.>'),end

% ----------------------------------------------------
% * set up the gui figure
% ----------------------------------------------------
f = figure('Color','w');
tmp = get(gcf,'Position') .* [1 1 .5 .1];
set(gcf,'Position',tmp)
set(gcf,'MenuBar','none','NumberTitle','off')
figure(f), set(gca,'Xlim',[0 100])
drawnow



% ----------------------------------------------------
% * create subdirectory
% ----------------------------------------------------
dirindex = 1;
if dirindex < 10, myz = '000';, else, myz = '00';, end
mydir = ['hewma' myz num2str(dirindex)];
eval(['mkdir ' mydir])
cd(mydir)




% ----------------------------------------------------
% * set up vars to save and save setup
% ----------------------------------------------------

% variables containing timeseries values
global grpmean
global grpt
global grpste
global grpp
global grph

% save setup stuff
try
    SETUP.ewma_stat_files = EXPT.FILES.ewma_z;
    SETUP.var_files = EXPT.FILES.ewma_var;
    SETUP.covariates = covt;
    SETUP.V = V;
    save SETUP SETUP
catch
    warning('Error creating SETUP file');
end

% ----------------------------------------------------
% * run robust regression for selected sets
% ----------------------------------------------------
    %if dools, disp('Running OLS and IRLS comparison (slower) - to turn this off, use 0 as 3rd input argument.'),end
    disp('____________________________________________________________')
    

for i = wh
    % * run robust HLM for each slice 
    % ----------------------------------------------------
    warning off
    fprintf(1,'Slice %3.0f. > ',i);

    if dools
        rob_fit(EXPT.FILES.ewma_z,EXPT.FILES.ewma_var,covt,i,dools,domask);
    else
        rob_fit(EXPT.FILES.ewma_z,EXPT.FILES.ewma_var,covt,i,dools,domask);
    end
    warning on
    fprintf(1,'\n');
end


cd ..



return



function [newP,newP2] = rob_fit(zimgs,varimgs,covt,index,dools,domask)

global xdim
global ydim
global emptyimage
global V
global domask

global grpmean
global grpt
global grpste
global grpp
global grph

global lam
global f

% --------------------------------------------
% load data
% --------------------------------------------
t1 = clock; fprintf(1,'Loading data. ');

L = size(zimgs{1},1);   % number of slices
if index > L, error('Slice index > number of slices.');,end
N = length(zimgs);      % number of subjects

for i = 1:N
    tmp = load(deblank(zimgs{i}(index,:))); % the slice ewma stat data
    z(:,:,i) = full(tmp.zdat);
    
    tmp = load(deblank(varimgs{i}(index,:))); % the slice variance data
    v(:,:,i) = full(tmp.vardat);
end
    
tp = size(z,2); % time points

fprintf(1,'%3.0f s. ',etime(clock,t1));
    
% --------------------------------------------
% find the in-analysis voxels
% --------------------------------------------

% mask, if spc
if domask
    Vm = spm_vol(domask);, vm = spm_read_vols(Vm);,

    % reshape mask for this slice
    vm = squeeze(vm(:,:,index));
    vm = reshape(vm,prod(size(vm)),1);
    vm = repmat(vm,[1 tp N]);
    
    z = z .* vm;
    v = v .* vm;
end

wh=sum(z,2);
wh = squeeze(wh);
wh = ~any(wh == 0 | isnan(wh),2);
wh = find(wh);      % which voxels are non-zero, non-nan for all subjects

tmp = length(wh); if tmp == 0, disp('No voxels in analysis!'), end


% --------------------------------------------
% set up output arrays
% --------------------------------------------

% to be written in .img files
tvals = NaN .* zeros(size(z,1),1);
pvals = NaN .* zeros(size(z,1),1);
sbmean = NaN .* zeros(size(z,1),1);
hvals = NaN .* zeros(size(z,1),1);
zvals = NaN .* zeros(size(z,1),1);
cp = NaN .* zeros(size(z,1),1);

grpmean{index} = zeros(size(z(:,:,1))) .* NaN;     % save betas
grpt{index} = zeros(size(z(:,:,1))) .* NaN;    % save t values
grpp{index} = zeros(size(z(:,:,1))) .* NaN;    % save p values
grpste{index} = zeros(size(z(:,:,1))) .* NaN;    % save ste values
grph{index}  = zeros(size(z(:,:,1))) .* NaN;    % save hyp test indicator
        

fprintf(1,'. %6.0f voxels. ',length(wh))


% --------------------------------------------
% perform regression
% --------------------------------------------
et = clock;

for i = 1:length(wh)
    vindx = wh(i);
    
    dat = squeeze(z(vindx,:,:))';
    vdat = squeeze(v(vindx,:,:))';

    [p,tm,Zcor,sbm,Zpop,ttime,sb,stats] = hewma(dat,vdat, lam);
    
    grpmean{index}(vindx,:) = Zpop;     % group timeseries
    grpt{index}(vindx,:) = ttime;       % t-values for each timepoint
    %grpp{index}(vindx,:) = p;
    grpste{index}(vindx,:) = sb;        % between-subjects standard dev (sigma)
    %grph{index}(vindx,1) = p<.05;       % global p-value
    
    % summary overall statistics
    tvals(vindx,1) = tm;
    pvals(vindx,1) = p;
    sbmean(vindx,1) = sbm;           % between-subjects st. dev (sigma), pooled over time
    hvals(vindx,1) = p < .05;
    zvals(vindx,1) = Zcor;
    tthresh(vindx,1) = stats.tthresh;
    cp(vindx,1) = stats.cp;
    
        %need corrected sig, change point, num total sig;
        %[m,t,p,se] = robust_mean(dat);
         
    if rem(i,10) == 0
        try,
            figure(f), 
            try,barh(100*i / length(wh)),catch,end
            set(gca,'Xlim',[0 100]),set(gca,'YTickLabel',i),drawnow
            text(5,.5,['Voxel ' num2str(i) ' of ' num2str(length(wh))],'Color','r')
        catch
        end
    end
   
end
fprintf(1,'\tDone in %3.0f s\n',etime(clock,et))

save hewma_timeseries grpmean grpt grpste lam xdim ydim

%reshape
tvals = reshape(tvals,xdim,ydim);
pvals = reshape(pvals,xdim,ydim);
sbmean = reshape(sbmean,xdim,ydim);
hvals = reshape(hvals,xdim,ydim);
zvals = reshape(zvals,xdim,ydim);
tthresh = reshape(tthresh,xdim,ydim);
cp = reshape(cp,xdim,ydim);

%write image slices

Pt = write_beta_slice(index,V,tvals,emptyimage,'hewma_t');
Pp = write_beta_slice(index,V,pvals,emptyimage,'hewma_p');
Ps = write_beta_slice(index,V,sbmean,emptyimage,'hewma_s');
Ph = write_beta_slice(index,V,hvals,emptyimage,'hewma_sig');
Pz = write_beta_slice(index,V,zvals,emptyimage,'hewma_z');
Pthr = write_beta_slice(index,V,tthresh,emptyimage,'hewma_thresh');
Pc = write_beta_slice(index,V,cp,emptyimage,'hewma_cp');


return






function Pw = write_beta_slice(slicei,V,betas,emptyimg,varargin)
% Pw = write_beta_slice(sliceindex,V,data,empty,prefix)
% Slice-a-metric version
warning off % due to NaN to int16 zero conversions
V.dim(4) = 16; % set to float to preserve decimals

prefix = 'image_';
if length(varargin) > 0, prefix = varargin{1};,end

for voli = 1:size(betas,3) % for each image/beta series point
    %if voli < 10, myz = '000';, elseif voli < 100, myz = '00';, else myz = '000';,end
    %V.fname = [prefix myz num2str(voli) '.img'];
    V.fname = [prefix '.img'];
    V.descrip = ['Hewma output image ' num2str(voli)];

    % create volume, if necessary
    if ~(exist(V.fname) == 2), spm_write_vol(V,emptyimg);,end
        
    spm_write_plane(V,betas(:,:,voli),slicei);
end

Pw = which(V.fname);

warning on
return

