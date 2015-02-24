function  PREPROC = canlab_task_independent_realignment(PREPROC, X, varargin)
% Augmented motion-correction that removes task-related activity before estimating realignment parameters for each image.
%
% Usage:
% -------------------------------------------------------------------------
% PREPROC = canlab_task_independent_realignment(PREPROC, X, [optional inputs])
%
% Author and copyright information:
% -------------------------------------------------------------------------
%     Copyright (C) 2013 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Inputs:
% -------------------------------------------------------------------------
% PREPROC       A struture created by canlab_preproc_2012 or canlab_preproc
% X             A design matrix for effects you want to remove before
%               calculating realignment parameters
% 'plot'        Create optional plots
% 'check'       Additional plots of pre- and post-realignment results
% 'noresid'     Do not residualize images; intended primarily as a check of the method
%
% Outputs:
% -------------------------------------------------------------------------
% xxx           A 4-D image of realigned data for the entire subject (multiple runs)
%
% Examples:
% -------------------------------------------------------------------------
%
% give examples here
%
% See also:
% canlab_preproc_2012

% Programmers' notes:
% List dates and changes here, and author of changes

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% Defaults
% -----------------------------------
% initalize optional variables to default values here.
%rowsz = [];
doplot = 0;
%basistype = 'spm+disp';
docheck = 0;
tinit = tic;
dosave = 0;
doresid = 1;  % residualize

% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            %case 'rows', rowsz = varargin{i+1}; varargin{i+1} = [];
            case 'check', docheck = 1;
            case 'plot', doplot = 1;
            case 'save', dosave = 1;
            
                case 'noresid', doresid = 0;
                    
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% Load data and ID outliers
% -------------------------------------------------------------------------

print_header1('CANlab Task-independent relignment');

% Load data for subject
dat = fmri_data(char(PREPROC.func_files{:}));

% Check X and Preproc
n = size(X, 1);
if n ~= size(dat.dat, 2)
    error('Number of images does not match number of rows in design (X) matrix!');
end

% View images we need:
% -----------------------------------
if doplot
    create_figure('X matrix to remove before realignment');
    plot_matrix_cols(X, 'horiz')
    canlab_preproc_montage_first_volumes(PREPROC.func_files);
end

% ID outliers (per session)
% -----------------------------------
% This duplicates what is done in canlab_preproc_2012
dat.images_per_session = PREPROC.images_per_session';
dat = preprocess(dat, 'outliers', 'plot');              % Spike detect and globals by slice
dat = preprocess(dat, 'outliers_rmssd', 'plot');        % RMSSD Spike detect

outliers = dat.covariates;

if docheck
    dat.Y = X;
    disp('X-matrix t-maps BEFORE TIR procedure.')
    [~, tmaps] = regress(dat, .001, 'unc', 'nodisplay');
end


% -------------------------------------------------------------------------
% PREPARE RESIDUALIZED DATASET FOR REALIGNMENT
% -------------------------------------------------------------------------

% Get path and base file name
outdirbase = fullfile(PREPROC.basedir, 'Functional', 'Preprocessed');


% Remove task data: temporary residual data for realign
% -----------------------------------

dat.covariates = X;  % this will be removed

if doresid
    resid_dat = preprocess(dat, 'resid', 1);  % flag: no add intercept
else
    resid_dat = dat;
end

resid_dat.fullpath = fullfile(outdirbase, 'tmp_resid_data.img');
write(resid_dat);

% -------------------------------------------------------------------------
% INTERPOLATE OUTLIERS
% -------------------------------------------------------------------------
% Optional.
% This will avoid spreading outlier points in time if slice timing is done
% later. Done after cov regression to avoid changing motion parameter
% estimates.

% But SHOULD use compressive sampling for this!!!!

% whole images can use interp1
% fancier methods should use interp3 or compressive sampling
%Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)

whout = sum(outliers, 2) > 0 ;
dat = preprocess(dat, 'interp_images', whout);

% -------------------------------------------------------------------------
% REALIGN AND WRITE
% -------------------------------------------------------------------------

% Realign
% -----------------------------------

flags = struct('quality',1,'fwhm',5,'sep',4,'interp',2,'wrap',[0 0 0],'rtm',0,'PW','','graphics',0,'lkp',1:6);

fprintf('Realigning. Reg to Mean = %3.0f. ', flags.rtm);
tic

% P will contain spm_vol vector with realigned mat estimates in .mat
P = spm_realign(resid_dat.fullpath, flags);
fprintf('Done in %3.0f sec\n', toc);

% Reslice (write realigned)
% -----------------------------------

% Use mat estimates from residual dataset in orignal (to realign) dataset
Pr = spm_vol(dat.fullpath);

for i = 1:length(Pr)
    
    Pr(i).mat = P(i).mat;  % use mats from realigned residual imgs
end

% check
if docheck
    datpre = fmri_data(dat.fullpath);
end

%ID images in each run
I = logical(intercept_model(PREPROC.images_per_session));

for i = 1:size(I, 2)  % run each run
    
    %wh = find(I(:, i));
    Prun = Pr(I(:, i));   % I(:, i) is which 3-D images
    
    save_parameters(Prun);
    
end

% reslice and write with new mats
% -----------------------------------

rflags = struct('interp',1,'mask',1,'mean',1,'which',2,'wrap',[0 0 0]', 'prefix','TIR_');
spm_reslice(Pr, rflags);


% Clean up
% -----------------------------------

% Remove temporary
torem = [resid_dat.fullpath(1:end-4) '.*'];
delete(torem)
% Move realigned files to Preprocessed
% Add names to PREPROC
PREPROC = canlab_preproc_clean_up_and_move_files(PREPROC);

% save realignment files and movement covariates

rpTIR = {};
nruns = length(PREPROC.TIR_mvmt_param_files);
for i = 1:nruns
    rpTIR{i} = importdata(PREPROC.TIR_mvmt_param_files{i});
end
[mvmtcat, mvmt_by_session, motion_cov_set] = canlab_preproc_motion_covariates(PREPROC.TIR_mvmt_param_files, PREPROC.images_per_session);

%generate temporary dat to adjust movement parameters
mdat.covariates = motion_cov_set;
mdat.images_per_session = PREPROC.images_per_session;
%split parameters into runs without removing any of them
mdat = split_and_clean_nuisance_covs(mdat,0);
%save output values to PREPROC
PREPROC.Nuisance.motion_cov_set = mdat.covariates;

if dosave
    disp('Saving results in updated PREPROC struture in PREPROC.Nuisance.motion_cov_set.')
    
    preprocmatname = fullfile(basedir, 'Functional', 'Preprocessed', 'PREPROC_SETUP.mat');
    savePREPROC(PREPROC, preprocmatname);
else
    disp('NOT saving PREPROC with updated realignment params. Enter ''save'' option to do so.');
end

% !rm realigned_data.img
% !rm realigned_data.hdr
% !rm realigned_data.mat
% !rm meanrealigned_data.img
% !rm meanrealigned_data.hdr

fprintf('REALIGN ALL DONE. Total time is %3.0f sec\n', toc(tinit));

% check
% -----------------------------------
if docheck
    print_header1('Checks: Compare standard realign with TIR');
    
    datpost = fmri_data(char(PREPROC.TIR_func_files{:}));
    
    create_figure('res', 2, 2);
    hist(mean(datpre.dat, 2), 100); title('PRE-means');
    subplot(2, 2, 2); hist(mean(datpost.dat, 2), 100); title('POST-means');
    subplot(2, 2, 3); hist(std(datpre.dat, [], 2), 100); title('PRE-std');
    subplot(2, 2, 4); hist(std(datpost.dat, [], 2), 100); title('POST-std');
    
    s1 = std(datpre.dat, [], 2);
    s2 = std(datpost.dat, [], 2);
    create_figure('Std'); plot(s1, s2, 'k.');
    xlabel('Pre'); ylabel('Post'); plot([0 250]', [0 250]', 'r--');
    title('Standard deviation across time at each voxel');
    
    m1 = mean(datpre.dat, 2);
    m2 = mean(datpost.dat, 2);
    create_figure('Means'); plot(m1, m2, 'k.');
    xlabel('Pre'); ylabel('Post'); plot([0 250]', [0 250]', 'r--');
    title('Mean at each voxel');
    
    try
        % may not work...
        s = mean(datpre);
        s.dat = [std(datpre.dat')' std(datpost.dat')']; % std across time
        s.dat(:, 3) = s.dat(:, 1) - s.dat(:, 2);        % difference in std
        s.dat = [s.dat mean(datpre.dat')' mean(datpost.dat')'];
        s.dat(:, 6) = s.dat(:, 4) - s.dat(:, 5);        % difference in mean
        orthviews(s)
        spm_orthviews_name_axis('PRE-std', 1);
        spm_orthviews_name_axis('POST-std', 2);
        spm_orthviews_name_axis('PRE-POST std', 3);
        spm_orthviews_name_axis('PRE-mean', 4);
        spm_orthviews_name_axis('POST-mean', 5);
        spm_orthviews_name_axis('PRE-POST mean', 6);
        
    catch
        disp('orthviews display failed....likely wrong sizes.');
    end
    
    % Viz movement params
    rp = {};
    nruns = length(PREPROC.mvmt_param_files);
    for i = 1:nruns
        rp{i} = importdata(PREPROC.mvmt_param_files{i});
    end
    
    if ~isempty(rp)
        canlab_preproc_motion_covariates(PREPROC.mvmt_param_files, PREPROC.images_per_session, 'plot');
        set(gcf, 'Tag', 'Movement params with standard realignment');
        
        try
            rpdiffs = cat(1, rp{:}) - cat(1, rpTIR{:});
            create_figure('Diffs: Standard realignment - TIR')
            plot(rpdiffs)
        catch
            disp('RP and RPTIR are different sizes!!!');
        end
        
        canlab_preproc_motion_covariates(PREPROC.TIR_mvmt_param_files, PREPROC.images_per_session, 'plot');
        set(gcf, 'Tag', 'Movement params with TIR');
        

        
%         [cc,stats] = cancor([X(:, 1:6) cat(1, rp{:})], 6);
%         [A,B,R,U,V,STATS] = canoncorr(X(:, 1:6), cat(1, rp{:}));
    end
    
end

end % function

% from SPM
function save_parameters(V)
fname = [spm_str_manip(prepend(V(1).fname,'rpTIR_'),'s') '.txt'];
n = length(V);
Q = zeros(n,6);
for j=1:n,
    qq     = spm_imatrix(V(j).mat/V(1).mat);
    Q(j,:) = qq(1:6);
end;
save(fname,'Q','-ascii');
end

function PO = prepend(PI,pre)
[pth,nm,xt,vr] = spm_fileparts(deblank(PI));
PO             = fullfile(pth,[pre nm xt vr]);
end


% get_basename(PREPROC.func_files{1}(1, :))
%
%     function get_basename(fn)
%
%         % strip ,num if exists
%         wh = find(fn == ',');
%         fn(wh:end) = [];
%
%         % get name, dir
%         [dd, basefile, ext] = fileparts(fn);
%         basefile = [basefile ext];
%
%     end

% From canlab_preproc_2012
function dat = split_and_clean_nuisance_covs(dat,doclean)

if nargin < 2
    %This eliminates doubled spikes.  The default is 'on'
    doclean = 1;
end
%create a copy of the covariates
spikes = dat.covariates;
%figure out how many images in each run
n_imgs = dat.images_per_session;
%get the number of runs
n_runs = numel(n_imgs);
%initialize the new covariate structure.  The first cells will be run
%specific and the final cell will contain covariates for use with
%concatenated data.
dat.covariates = cell(n_runs+1,1);
for i = 1:n_runs
    %choose the rows that contain data for this run and get those covs
    wh_rows = (1:n_imgs(i))+sum(n_imgs(1:i-1));
    covs = spikes(wh_rows,:);
    %if cleaning is requested,
    %check for copies in these spike regressors.
    %This may occur since spikes are defined both by Mahalanobis
    %distance and standard deviation
    covs = covs(:,var(covs)~=0);
    n_covs = size(covs,2);
    if doclean
        cov_pos = zeros(n_covs,1);
        for j = 1:n_covs
            cov_pos(j,1) = find(covs(:,j));
        end
        %only keep unique covs
        [dummy wh_covs] = unique(cov_pos);
    else
        wh_covs = 1:n_covs;
    end
    dat.covariates{i} = covs(:,wh_covs);
end

%repeat above loop, but use the entire concatenated set of covariates
n_spikes = size(spikes,2);
if doclean
    spike_pos = zeros(n_spikes,1);
    for i=1:n_spikes
        spike_pos(i) = find(spikes(:,i));
    end
    [dummy wh_spikes] = unique(spike_pos);
else
    wh_spikes = 1:n_spikes;
end
dat.covariates{n_runs+1} = spikes(:,wh_spikes);

%do a quick check to make sure that all the rows sum together properly
all_scans = size(dat.covariates{n_runs+1},1);
run_scans = 0;
for i = 1:n_runs
    run_scans = run_scans + size(dat.covariates{i},1);
end
if all_scans ~= run_scans
    error('The nuisance covariate processing has encountered an error: The number of scans per run do not match.')
end

end


function savePREPROC(PREPROC, preprocmatname)
fprintf('Saving PREPROC setup structure with image names, etc., in file: \n%s\n', preprocmatname)

if exist(preprocmatname, 'file')
    try
        save(preprocmatname, '-append', 'PREPROC')
    catch
        % may not be able to find; perhaps doesn't exist, but file does?
        disp('Warning: Overwriting PREPROC because PREPROC variable not found in .mat file.')
        save(preprocmatname, 'PREPROC')
    end
else
    save(preprocmatname, 'PREPROC')
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

