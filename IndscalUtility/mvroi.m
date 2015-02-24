 function DATA=mvroi(varargin)
% function DATA=mvroi([DATA],[SPEC])
%
% INPUTS:
%
% DATA can be a structure with DATA.dat, dat, or the (string) name of a file containing dat
%   dat     data, cell vector, one cell per subject, time x regions matrix
%           in each cell
%
% SPEC can be a structure with options for the analysis or the name of a
%           file containing the SPEC structure (called SPEC)
%           SPEC must have at least event onsets in SPEC.onsets or model matrices for each
%           subject in SPEC.DX, or no task will be removed
%
%
% MVROI - a multivariate analysis tool for miniblock fMRI designs 
% type help mvroi_specs for more information
%
%
% Examples:
% If you have data (dat) saved in input_data.mat, and specs (SPEC) saved in
% mvroi_specs.mat (recommended name):
% DATA = mvroi('input_data','mvroi_specs');
%
% If you have an existing DATA structure with some steps completed:
% load DATA
% DATA = mvroi(DATA,'mvroi_specs');
%
close all



diary mvroi_output.txt



% -----------------------------------------------------------------------
% STEP 1 - load data,specs & names
% -----------------------------------------------------------------------
% Checks whether a spec file exists in the path (prefix_specs.mat)
% If not, generates a spec file through inputs.
% calls mvroi_specs

if length(varargin)==0
    help mvroi
    error('No input arguments!')
    
else
    % define DATA
    tmp = varargin{1};
    if isstr(tmp), 
        load(tmp);,
        try, DATA.DATA.dat = dat;, catch, error('Your file must contain a variable called dat.'), end
    elseif iscell(tmp)
        DATA.DATA.dat = tmp;
    elseif isstruct(tmp)
        DATA = tmp;
        if ~isfield(DATA.DATA,'dat'), error('Your DATA structure must contain a field called DATA.DATA.dat.'),end
    end
end

if length(varargin) > 1
    % define SPEC
    tmp = varargin{2};
    if isstr(tmp), 
        load(tmp);,
        try, SPEC = SPEC;, catch, error('Your 2nd arg file must contain a variable called SPEC.'), end
    elseif isstruct(tmp)
        SPEC = tmp;
    end
else
    
    % we must define a new SPEC structure
    SPEC = [];

end

% Check existing fields in SPEC, and prompt to enter any that are missing
DATA.SPEC = mvroi_specs(SPEC,DATA);



disp('Saving DATA.mat'), save DATA DATA





% -----------------------------------------------------------------------
% STEP 2 - generate design matrices
% -----------------------------------------------------------------------
% (If they don't aready exist)
% Check for event onsets
% Replicate first subject's onsets and design matrix if only one exists
% (Assume same design for each subject in this case)
%
% Needs:
% DATA.SPEC.onsets      onsets of event times in SPM-style format (cell array
%                       of times, condition nested within session
% numsub                number of subjects
% DATA.SPEC.firpoints   vector or integer of number of betas to estimate
%                       for each condition; same as EXPT.DX.numframes
% DATA.SPEC.spersess    Number of images in each session, [n1 n2 ... nx]
%
%
% Will take:
% DATA.DX        supersedes DATA.onsets
%
% Creates:
% DATA.DX        cell array of one design matrix per subject
%
% Uses functions (nested to reflect structure of calls):
% build_fir_model
%   onsets2delta
%   tor_make_deconv_mtx3
%       intercept_model
%       

disp(' ')
disp('% -----------------------------------------------------------------------')
disp(' STEP 2 - Generating design matrices for each subject')
disp('% -----------------------------------------------------------------------')
mvroi_generate_design_plugin

if dobuild, disp('Saving DATA.mat'), save DATA DATA, end




% -----------------------------------------------------------------------
% STEP 3 - filter data, remove task-related variance
% -----------------------------------------------------------------------

%%%-------scaling, high-pass, trimming----------%%%
% experiment-critical fields for filtering, with example values:
% DATA.SPEC.trim = 3, number of standard deviations to trim to
% DATA.SPEC.HP = 120;
% DATA.SPEC.TR = 1.5;
% % output is DATA.DATA.filtered_dat

disp(' ')
disp('% -----------------------------------------------------------------------')
disp(' STEP 3 - filter data, remove task-related variance')
disp('% -----------------------------------------------------------------------')


if isfield(DATA.DATA,'filtered_dat')
    dofilt = input('DATA.DATA.filtered_dat exists.  Re-trim and high-pass filter (1/0)? ');
else
    dofilt = 1;
end

if dofilt
    fprintf(1,'Preprocessing . TR = %3.2f, trim = %3.0f, HP = %3.0f\n',DATA.SPEC.TR,DATA.SPEC.trim,DATA.SPEC.HP);
    DATA=scaledata(DATA);
    dat=DATA.DATA.dat;
    disp('Number of points windsorized, Subjects x Regions:');
    print_matrix(DATA.DATA.ntrimmed,SPEC.names);
end

%%%---------------------------------------------%%%


%%%---- calculate betas, fit and residuals -----%%%
% Inputs:
% Requires DATA.DATA.filtered_dat field
%
% Outputs:
% DATA.DATA.resids, residuals of filtered data
% DATA.DATA.fits, timeseries of model fits for each region
% DATA.DATA.b, vector of FIR estimates concatenated
% DATA.DATA.fir is FIR HRF estimates for each event type; cells are {regions x event types},
% elements are (subjects x time)

if isfield(DATA.DATA,'resids')
    domodel = input('DATA.DATA.resids exists.  Re-remove model from filtered data (1/0)? ');
else
    domodel = 1;
end

if domodel
    DATA=mvroi_fit_model(DATA);
    
    try, saveas(gcf,'step3_hrf_estimates','fig'),saveas(gcf,'step3_hrf_estimates','tif'), close, catch, disp('Trouble saving hrf_estimates figure.'),end 
end


%%%---------------------------------------------%%%


%%%--remove global timeseries by regression ----%%%
%
% input:  DATA.resids
% output: DATA.resids (global mean removed)
%   also returns DATA.DATA.globalmeans
%
% Only do this if we're removing task, so we don't do it twice if we skip
% model removal stage

if domodel
    
    switch DATA.SPEC.doglobal
        case 'none'
            DATA.DATA.resids_descrip = 'No global scaling or mean adjustment';
        case 'regression'
            DATA=mvroi_remove_global_mean(DATA);
        otherwise
            disp('Warning!!  DATA.SPEC.doglobal should be ''regression'' or ''none'' . No global scaling/regression done');
            DATA.DATA.resids_descrip = 'No global scaling or mean adjustment';
    end
end
%%%---------------------------------------------%%%


%%%----plot hrfs and residuals by region --------%%%
doresplot = 0;
if doresplot,
    mvroi_residplot_plugin
end
%%%---------------------------------------------%%%

if domodel, disp('Saving DATA.mat'), save DATA DATA, end





% -----------------------------------------------------------------------
% STEP 3.5 - trial_level betas
% -----------------------------------------------------------------------
disp(' ')
disp('% -----------------------------------------------------------------------')
disp(' STEP 3.5 - Trial-level betas -- in development!')
disp('% -----------------------------------------------------------------------')
mvroi_trial_beta_plugin







% -----------------------------------------------------------------------
% STEP 4 - find number of non-gaussian dimensions & reduce dimensions
% -----------------------------------------------------------------------

disp(' ')
disp('% -----------------------------------------------------------------------')
disp(' STEP 4 - Choose number of dimensions')
disp('% -----------------------------------------------------------------------')


%%%---------decide whether to run choose_ndims-----------%%%
if strcmp(char(DATA.SPEC.ndims),'choose')  %don't run this if we hardcoded the dimensionality  
    dochoose = 1;
else
     %%%------- we already have some -------------%%%   
     if isstr(DATA.SPEC.ndims), DATA.SPEC.ndims = str2num(DATA.SPEC.ndims);,end
     dochoose = input(['Found ' num2str(DATA.SPEC.ndims) ' dimensions in DATA.SPEC.ndims.  Re-run (1/0)?']);   
end    
    

if dochoose
    %%%---------pick number of dimensions-----------%%%
        tor_fig;
    
        [ndims,DATA.NDIMS] = choose_ndims(DATA.DATA.resids);     % nonpar thresholding for significant pcs - on cov mtx
        if ndims == 0, ndims = 2;, warning('NO significant eigenvalues!  May be independent regions.  Expanding to 2-D space for display anyway.'),end
        if ndims == 1, ndims = 2;, warning('1-D Solution!  Expanding to 2-D space anyway.'),end
        fprintf(1,'\nChosen %3.0f dimensions.\n',ndims)
        DATA.SPEC.ndims = ndims;
        
        %DATA.SPEC.ndims = input('Enter number of dims to use: ');%%ask anyway!
        
        try, saveas(gcf,'step4_eigenvalues','tif');,
            saveas(gcf,'step4_eigenvalues','fig');
        catch disp('cannot save figure');,
        end
        
        % save SPEC structure
        disp('Saving mvroi_specs.mat'),SPEC = DATA.SPEC; save mvroi_specs SPEC
        
        disp('Saving DATA.mat'), save DATA DATA

    %%%---------------------------------------------%%%           
end













% -----------------------------------------------------------------------
% STEP 5 - get cross-lagged correlations on residuals
% -----------------------------------------------------------------------
%%%% correlate betas, fit and residuals %%%%%%%%
% Input: DATA.DATA.resids
% Outputs:
%   DATA.DATA.corr     cell vector of correlation matrices for each subject
%   DATA.DATA.lat      cell vector of latency matrices for each subject
%   DATA.DATA.xc       3-D matrix, reg x reg x state within subject
%                       suitable for input to decomposition 

% we did the regressoutmean above, don't do it again.  DO it above so it's
% done before choosing number of dimensions.

disp(' ')
disp('% -----------------------------------------------------------------------')
disp(' STEP 5 - cross-correlations')
disp('% -----------------------------------------------------------------------')

if isfield(DATA.DATA,'xc')
    doxc = input('DATA.DATA.xc exists.  Re-do cross correlations (1/0)? ');
else
    doxc = 1;
end

if doxc
    [DATA,xc,xl,resid] = mvroi_xc(DATA,'doresid','no','doglobal','none');
    %[DATA,xc,xl,resid] = mvroi_xc(DATA,'doresid','no','doglobal','regression');

    %DATA.meanxc = mean(DATA.xc,3);  % mean across tasks and subjects
    %DATA.meanxl = mean(DATA.xl,3);

    disp('Saving DATA.mat'), save DATA DATA
end













% -----------------------------------------------------------------------
% STEP 6 - multivariate data decomposition
% -----------------------------------------------------------------------
%%  run indscal & plot initial group solution %%%%%%%
%everything from this step and from the clustering is saved in DATA.INDSCAL
%
% indscal3 returns Group Space (Gs) of stim coordinates, Weight matrix W,
% a stats structure (see help in indscal3) and scalar product matrices.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% JUST IN CASE WE FORGET                           %%
%% a DISSIMILARITY matrix has zeros on the diagonal %%
%% a SIMILARITY matrix has ones on the diagonal     %%
%% neither can have negative values                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('% -----------------------------------------------------------------------')
disp(' STEP 6 - INDSCAL')
disp('% -----------------------------------------------------------------------')

if DATA.SPEC.betaflag
    % do on beta weights -- we have to scale them so we get valid sim
    % matrix
    xc = DATA.DATA.xc; xc = xc ./ max(abs(xc(:)));  % scale to max of 1; arbitrary...
    xc = (1 - xc) ./ 2;
    
    % make sure zeros are on diagonals
    mmat = 1 - eye(size(xc,1),size(xc,2));
    mmat = repmat(mmat,[1 1 size(xc,3)]);
    xc = xc .* mmat;
    DATA.DATA.dissim = xc;
    
else
    % scale correlations (or beta weights) so that they are a similarity matrix between 0 and 1.
    DATA.DATA.dissim = (1 - DATA.DATA.xc) ./ 2;
end

%[Gs,W]=indscal1(DATA.rdat,ndims,names);
[Gs,W,DATA.INDSCAL.indscalstats,DATA.DATA.scalar_products]=indscal3(DATA.DATA.dissim,DATA.SPEC.ndims,DATA.SPEC.names,DATA.SPEC.corrscaling);

DATA.INDSCAL.Gs = Gs;
DATA.INDSCAL.W = W;

disp('Saving DATA.mat'), save DATA DATA

% stats and plots
mvroi_weight_test_plugin

disp('Saving DATA.mat'), save DATA DATA

DATA.INDSCAL.indiv_weights = mvroi_individual_weights(DATA);






% -----------------------------------------------------------------------
% STEP 7 - clustering
% -----------------------------------------------------------------------

disp(' ')
disp('% -----------------------------------------------------------------------')
disp(' STEP 7 - clustering of regions')
disp('% -----------------------------------------------------------------------')


docluster = 1; 
if isfield(DATA,'CLUSTER'), docluster = input('Existing CLUSTERS found.  Recompute cluster assignments (1/0)? ');,end

if docluster

numreg = DATA.SPEC.numreg;  %number of ROIs
    
%%%---- get number of clusters to search over -----%%%
if isstr(DATA.SPEC.clustsolutions) 
    % if it's 'choose', converts to integers; otherwise, assumes it's already integers
   
    if ~strcmp(DATA.SPEC.clustsolutions,'choose')
        % assume it's numbers in string format
        DATA.SPEC.clustsolutions = str2num(DATA.SPEC.clustsolutions);
    else
        % it's 'choose'
        % automatically choose number of clusters to search over, if 'choose'
        % entered in SPEC
        DATA.SPEC.clustsolutions = 2:round((numreg - 1) ./ 2);
    end
    
end
fprintf(1,'Searching over k clusters: %3.0f to %3.0f . ',min(DATA.SPEC.clustsolutions), max(DATA.SPEC.clustsolutions))
%%%---------------------------------------------%%%


%%%---- run clustering -------------------------%%%
[DATA.CLUSTER.pval,DATA.CLUSTER.classes,DATA.CLUSTER.names,DATA.CLUSTER.Gs, ...
    DATA.CLUSTER.include,clustnames,DATA.CLUSTER.stats]= ...
    testclustnew(Gs,DATA.SPEC.clustsolutions,DATA.SPEC.ndims, ...
    DATA.SPEC.numperm,DATA.SPEC.names,DATA.SPEC.remove,'average');

% remove entries for nonselected regions;
DATA.CLUSTER.xcsave = DATA.DATA.xc;
DATA.CLUSTER.xcsave=DATA.CLUSTER.xcsave(DATA.CLUSTER.include,DATA.CLUSTER.include,:);
DATA.CLUSTER.meanxc=mean(DATA.CLUSTER.xcsave,3);
%%%---------------------------------------------%%%

disp('Saving DATA.mat'), save DATA DATA

try, saveas(gcf,'step7_testclust_details','tif');, saveas(gcf,'step7_testclust_details','fig');,catch disp('cannot save figure');,end

% plot group solution
if size(DATA.CLUSTER.Gs,2) <= 2, tor_fig;, end      % if 2D, make a figure, otherwise unnecessary
mdsfig(DATA.CLUSTER.Gs,DATA.CLUSTER.names,DATA.CLUSTER.classes);
drawnow
if size(DATA.CLUSTER.Gs,2) <= 2,
    try saveas(gcf,'step7_groupspace','tif');,catch disp('cannot save figure');,end
else
    % 3 or more dims
    try saveas(gcf,'step7_groupspace','tif');,close, 
        saveas(gcf,'step7_groupspace3d','tif');,close,
    catch disp('cannot save figure');,
    end
end


end % if docluster









% -----------------------------------------------------------------------
% STEP 8 - Stats on correlations (all regions)
% -----------------------------------------------------------------------
% print correlations, average across conditions and diff among conditions
% uses contrast3d and xcstats to print statistics for average and contrast
% All correlations for all regions, not just included ones.

disp(' ')
disp('% -----------------------------------------------------------------------')
disp(' STEP 8 - Stats on individual regional correlations')
disp('% -----------------------------------------------------------------------')

docorrel = 1; 
if isfield(DATA,'CORRELS'), input('Existing CORRELS found.  Recompute stats, output and plots (1/0)? ');,end

if docorrel

    DATA.CORRELS = mvroi_correl_stats(DATA,DATA.SPEC.names);

    % plot group solution and differences among conditions 

    mvroi_mdsfig_plot2(DATA.CLUSTER,DATA.SPEC,DATA.CORRELS.AVGSTATS.sigmat_uncorrected,DATA.CORRELS.DIFSTATS.sigmat_uncorrected,'Uncorrected');
 
    try, saveas(gcf,'step8_groupspace_lines_uncorrected','tif');,
        saveas(gcf,'step8_groupspace_lines_uncorrected','fig');
    catch disp('cannot save figure');,
    end

    mvroi_mdsfig_plot2(DATA.CLUSTER,DATA.SPEC,DATA.CORRELS.AVGSTATS.sigmat,DATA.CORRELS.DIFSTATS.sigmat,'Corrected');
 
    try, saveas(gcf,'step8_groupspace_lines_corrected','tif');,
        saveas(gcf,'step8_groupspace_lines_corrected','fig');
    catch disp('cannot save figure');,
    end

    disp('Saving DATA.mat'), save DATA DATA
    
    % plot separate states
    mvroi_mdsfig_plot_sepstates(DATA.CLUSTER,DATA.CORRELS,DATA.SPEC);
        try, saveas(gcf,'step8_sep_states_corrected','tif');,
        saveas(gcf,'step8_sep_states_corrected','fig');
        catch disp('cannot save figure');,
        end
    
end











% -----------------------------------------------------------------------
% STEP 9 - Name Networks, Stats, and Plot 
% -----------------------------------------------------------------------

% now apply the cluster solution
% DATA.SPEC.comps(1,:) is applied here as well.  
% as are stats on correlations between group/contrast values and behavioral
% regressors, if DATA.SPEC.beh is entered.




% first, name each network
name_networks_plugin

if ~isfield(DATA.SPEC,'beh'), DATA.SPEC.beh = [];,end

% get group correlations, lumping across task states
grpweights = ones(size(DATA.SPEC.comps(1,:)))./size(DATA.SPEC.comps,2); % to get average across conditions in apply_cluster_solution

DATA.APPLY_CLUSTER.group = apply_cluster_solution(DATA.CLUSTER.classes,DATA.DATA.xc,'state',[], ...
    'contrast',grpweights,'names',DATA.APPLY_CLUSTER.names,'bcov',DATA.SPEC.beh);

title('group correlations','FontSize',16);
try, saveas(gcf,'step9_group_correl_bar','tif');,saveas(gcf,'step9_group_correl_bar','fig');, close,     catch, disp('Error saving figure.'), end
try, saveas(gcf,'step9_group_correl_matrix','fig');,saveas(gcf,'step9_group_correl_matrix','tif'); close,     catch, disp('Error saving figure.'), end

if DATA.SPEC.numstates > 1
% get differences in correlations between task states
DATA.APPLY_CLUSTER.contrast = apply_cluster_solution(DATA.CLUSTER.classes,DATA.DATA.xc,'state',[], ...
    'contrast',DATA.SPEC.comps(1,:),'names',DATA.APPLY_CLUSTER.names,'bcov',DATA.SPEC.beh);

title(['contrast ',num2str(DATA.SPEC.comps(1,:))],'FontSize',16);
try, saveas(gcf,'step9_contrast_correl_bar','tif');, saveas(gcf,'step9_contrast_correl_bar','fig');, close,     catch, disp('Error saving figure.'), end
try, saveas(gcf,'step9_contrast_correl_matrix','fig');,saveas(gcf,'step9_contrast_correl_matrix','tif'); close,     catch, disp('Error saving figure.'), end

end % contrast


% get class (network) average stimulus locations (group space)
Gs = DATA.CLUSTER.Gs;

for i = 1:max(DATA.CLUSTER.classes)
    wh = find(DATA.CLUSTER.classes == i);
    if length(wh)>1;
    classGs(i,:) = mean(Gs(wh,:));
    else
    classGs(i,:) = Gs(wh,:);
    end
end
DATA.APPLY_CLUSTER.Gs = classGs;
DATA.APPLY_CLUSTER.classes = 1:max(DATA.CLUSTER.classes);


% Figures
mvroi_mdsfig_plot2(DATA.APPLY_CLUSTER,DATA.SPEC,DATA.APPLY_CLUSTER.group.group_stats.sigu,DATA.APPLY_CLUSTER.contrast.group_stats.sigu,'Uncorrected');
try, saveas(gcf,'step9_clusterspace_lines_uncorrected','tif');, saveas(gcf,'step9_clusterspace_lines_uncorrected','fig');, close,     catch, disp('Error saving figure.'), end

mvroi_mdsfig_plot2(DATA.APPLY_CLUSTER,DATA.SPEC,DATA.APPLY_CLUSTER.group.group_stats.sig,DATA.APPLY_CLUSTER.contrast.group_stats.sig,'Corrected');
try, saveas(gcf,'step9_clusterspace_lines_corrected','tif');, saveas(gcf,'step9_clusterspace_lines_corrected','fig');, close,     catch, disp('Error saving figure.'), end


% Tables

disp('Corrected average correlations')
m = DATA.APPLY_CLUSTER.group.mean_corrs; sig = abs(DATA.APPLY_CLUSTER.group.group_stats.sig);
n = DATA.APPLY_CLUSTER.names;
correlation_to_text(m,sig,n);

disp('Uncorrected average correlations')
m = DATA.APPLY_CLUSTER.group.mean_corrs; sig = abs(DATA.APPLY_CLUSTER.group.group_stats.sigu);
correlation_to_text(m,sig,n);

if DATA.SPEC.numstates > 1
    disp(['Corrected contrast across states, ' DATA.SPEC.comptitle])
    m = DATA.APPLY_CLUSTER.contrast.mean_corrs; sig = abs(DATA.APPLY_CLUSTER.contrast.group_stats.sig);
    correlation_to_text(m,sig,n);

    disp('Uncorrected contrast across states')
    m = DATA.APPLY_CLUSTER.contrast.mean_corrs; sig = abs(DATA.APPLY_CLUSTER.contrast.group_stats.sigu);
    correlation_to_text(m,sig,n);

end % contrast

disp('Saving DATA.mat'), save DATA DATA



close all;
diary off


return




