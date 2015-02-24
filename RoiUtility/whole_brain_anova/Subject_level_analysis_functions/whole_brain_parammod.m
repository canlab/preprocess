function Pw = whole_brain_parammod(P,DX,TR,HP,nruns,whpredictors,contrasts,connames,varargin)
    % function Pw = whole_brain_parammod(P,DX,TR,HP,nruns,whpredictors,contrasts,connames,[doplot],[mask],[smoothlen],[startslice])
    %
    % single-subject timeseries extraction and model fitting
    % saves filtered images f* and full-model betas (model fits)
    %
    % P = list of image names, raw functionals or beta images (see "type"
    %       below)
    % DX = model matrix with parametric modulators added to front, see
    % c2rtpredictor.m, intext_make_RT_models.m, rt2delta.m
    %
    % TR = repetition time of your study (sampling rate)
    % HP = high-pass filter cutoff in seconds to apply (SPM style)
    % nruns = number of sessions (intercepts) OR num. of images in each sess,
    %           e.g., [169 169 172]
    % whpredictors = indices of which predictors to write (save) beta images for
    %       could be parametric modulators, or param mods + events of interest
    % contrasts = matrix with one row per contrast.  should have same number of
    % columns as whpredictors.
    %
    % [doplot] = optional graphics
    % [mask] = optional name of mask image to limit processing to brain, or []
    % [smoothlen] = exponential smoothing filter length for timeseries and
    %           betas; influence is 0 after smoothlen seconds; default = 6
    % [startslice] = slice to start at, default = 1
    % [type] = optional type of analysis to run.  Options:
    %        'full'    Everything, from filtering to height/time/width, enter
    %                  raw image names and model matrix DX
    %
    %        'data'    Filter data only and save images in f*img
    %                  Enter raw image names and empty model matrix []
    %
    %        'htw'     Smooths beta series and computes height/time/width only
    %                  Enter beta images from FIR model within condition,
    %                   e.g., beta1_event1 beta2_event1 ... beta1_event2 b2_e2
    %                   ...
    %                  One possibility is to use images from SPM2 at this
    %                  stage.
    %
    % For a shell to run this function, see whole_brain_fir.m
    %
    % Tor Wager, 2/20/05
    %
    % %
    % The analysis sequence:
    % Run whole_brain_fir, which mainly just calls whole_brain_filter
    % whole_brain_filter gives you:
    %       1) trimmed, filtered images, saved in f*img
    %           these can be passed to SPM2, e.g., if desired
    %       2) extraction of beta images for FIR model, saved in dx_beta*img
    %           in individual subject folders; also smooths betas if desired
    %       3) calculation of height, delay, and width for each FIR extracted
    %           saved on images in individual subject folders.
    %
    % After running this, you'll need to collect image names for use in random
    % effects analyses, and most likely create contrasts across conditions for
    % differences in height, delay, and width.
    % to do this, first go to main model directory, with subject subfolders
    % to save a list of height/delay/width images:
    %EXPT = get_htw_image_names(EXPT)
    % to create contrasts across those images and save in EXPT.SNPM for rfx
    % analysis:
    %EXPT = make_htw_contrast_images(EXPT);
    % Now you can run robust regression across the images in EXPT.SNPM.
    % The robust reg uses image names in EXPT.SNPM.P, so you may have to save
    % images from another field (e.g., EXPT.SNPM.heightP) in the .P field.
    % you may also want to test individual contrasts against zero, in which
    % case you can use EXPT.NLCON.height images, for example (and others in
    % NLCON) instead.
    % EXPT = robfit(EXPT);


    % Explanations of some more variables
    % vb    = [optional] verbose output level: 0 none, 1 some, 2 lots
    % mask  = [optional] mask 3-D volume to apply before extracting
    %
    % P     = image file names (str matrix)
    % S     = filtering matrix (e.g., high-pass)
    % DX    = full model to fit, unfiltered
    % vb    = [optional] verbose output level: 0 none, 1 some, 2 lots
    % mask  = [optional] mask 3-D volume to apply before extracting
    % nsess = number of sessions (intercepts): assumes last nsess columns
    %         of DX are run intercepts to be removed before trimming
    %
    % dims  = dimensions of image files in data
    % cols  = no. of columns of DX, also no. of beta images to write
    % SDX   = smoothed (filtered) full model to fit
    % PSDX  = pseudoinverse of filtered full model
    % PSDXS = PSDX * S, ready to multiply with data y for pinv(SX) * Sy
    % betas = 4-D array of x,y,z,column
    %
    % ntrimmed = number of outliers trimmed from timeseries at each voxel
    % Pw    = string matrix of output file names
    %
    % Pw = whole_brain_filter(d,c.model(1:166,:),2,120,5,c.numframes,1,10);

    % -------------------------------------------------------------------
    % * configure input arguments
    % -------------------------------------------------------------------

    if isempty(P), warning('Empty image names!  No processing done.'), Pw = [], return, end

    type = 'full';    % predefined set of functions to do -- 'data' or 'full'
    if length(varargin) > 4, type = varargin{5};, end

    % defaults
    % ---------------
    % for data filtering
    global dofilt
    global writefilt

    % for extracting and parameterizing hrfs
    global doysmooth
    global dosmooth
    global dowritebetas
    global dohtw

    % for plotting
    global f1
    global f2
    global f3

    % for saving betas and contrasts
    global whpredictors
    global contrasts
    global connames

    % filtering
    % ---------------
    dofilt = 1;     % filter data
    writefilt = 0;  % write filtered data image set f*

    % model fitting
    % ---------------
    dohtw = 0;      % extract height, delay, and width
    dosmooth = 0;   % smooth beta series before hdw -- never do these things for canonical param mod model
    doysmooth = 0;  % smooth data before model fitting
    dowritebetas = 1;  % write beta images to disk as pmbeta_*img

    switch type
        case 'full', % do nothing
        case 'data', % data filtering only
            dohtw = 0; dosmooth = 0; doysmooth = 0; DX = [];writefilt = 1;
            fprintf(1,'\nWriting filtered data f*img only in:\n  %s\n',P(1,:))
        case 'htw', % enter beta images, do dosmooth and height/time/width only
            dofilt = 0; writefilt = 0;  doysmooth = 0;
            dohtw = 1; dosmooth = 1; dowritebetas = 0; DX = [];
        otherwise
            error('Unknown instruction type')
    end




    % high-pass smoothing
    smoothlen = 6;      % smoothing for beta series, in seconds (time to no effect)
    if length(varargin) > 2, smoothlen = varargin{3};, end


    % first slice -- option to start at something besides 1
    firstslice = 1;
    if length(varargin) > 3, firstslice = varargin{4};, end

    % plotting
    global dographics
    dographics = 1;
    if length(varargin) > 0, dographics = varargin{1};, end
    if dographics,  drawnow; end

    % verbose output
    vb = 2; %if length(varargin) > 0, vb = varargin{1};, end
    if vb > 0, t1 = clock;, fprintf(1,'\n\t\tSubject setup (mapping images): %s',pwd),end

    % get rid of empty image names
    wh = find(all(P == ' ',2));
    if ~isempty(wh),disp('WARNING: Some empty rows in image names.  Removing these.');
        P(find(wh),:) = [];
    end

    % image dimensions, etc.
    V = spm_vol(P(1,:));
    Vall = spm_vol(P);

    global dims
    global cols

    dims = V.dim(1:3);
    cols = size(DX,2);

    % initialize beta matrix
    betas = NaN * zeros([dims cols]);

    if vb > 0, fprintf(1,'\n\t\tNew image dims: %3.0f %3.0f %3.0f %3.0f ',dims(1), dims(2), dims(3), cols),end

    % brain masking
    mask = ones(dims);
    if length(varargin) > 1,
        mask = varargin{2};,
        if ischar(mask), Vm = spm_vol(mask);, mask = spm_read_vols(Vm);,end
    end
    if isempty(mask), mask = ones(dims);,end

    ntrimmed = NaN * zeros(dims);

    %if dographics, close, end

    % -------------------------------------------------------------------
    % * set filter
    % -------------------------------------------------------------------
    npoints = size(P,1);

    if dofilt

        fprintf(1,'\n\t\tSetting up filter for %3.0f images in %3.0f runs',npoints,length(nruns))

        % Get Residual-forming matrix hat such that hat * y  removes HP filter and
        % intercepts for each session; save in S

        [dummy,KL,KH] = use_spm_filter(TR,npoints,'none','specify',HP);		% HP filter mtx
        clear dummy KL
        X = intercept_model(nruns); % nruns has number of images in each run

        if size(X,1) ~= size(KH,1), error('Filter and intercept len do not match.  Wrong number of runs?'),end

        S = [KH X];                             % filtering matrix
        PS = pinv(S);
        hat = S * PS;                           % residual forming matrix
        S = eye(length(hat)) - hat;             % DX-hat*DX to get residuals = (I-hat)DX
        clear hat PS X

    else
        % no filtering
        fprintf(1,'\n\t\tNo data filtering, %3.0f images in %3.0f runs',npoints,nruns)
        S = eye(size(DX,1));

    end  % if dofilt

    if doysmooth
        [tmp,F] = smooth_timeseries(ones(npoints,1),smoothlen); % F is LP filter matrix for data
    else,
        F = [];
    end

    % -------------------------------------------------------------------
    % * adjust model DX based on filter, if not empty
    % -------------------------------------------------------------------

    if ~isempty(DX)
        fprintf(1,'...adjusting model matrix\n ')
        if size(S,2) ~= size(DX,1), error('Images and model size do not match!');, end

        SDX = S * DX;                           % model with intercept, HP terms removed, LP smooth
        PSDX = pinv(SDX);                       % pseudoinverse of filtered design matrix
        PSDXS = PSDX * S;                       % beta-forming matrix, post-multiplied by filter
        % so that PSDXS * y = betas
        if dographics,
            figure('Color','w'); subplot(1,2,1); set(gca,'FontSize',16)
            imagesc(S),title('Filtering matrix'), subplot(1,2,2); imagesc(SDX),title('Design matrix'),colormap hot, drawnow,
            saveas(gcf,'scnlab_filters','tif')
        end

    else
        SDX = []; PSDX = []; PSDXS = [];
    end

    if vb > 0, fprintf(1,'\n\t\tFinished in %3.0f s\n',etime(clock,t1)),end


    % -------------------------------------------------------------------
    % * for each slice...
    % -------------------------------------------------------------------

    for slicei = firstslice:dims(3)

        if vb > 0, t1 = clock;, fprintf(1,'\nSlice %3.0f \n------------>\n ',slicei),end

        [a,b] = process_slice(slicei,P,S,PSDXS,TR,V,Vall,F,mask(:,:,slicei),smoothlen);

        ntrimmed(:,:,slicei) = b;

        %dographics = 0;
        if vb > 0, fprintf(1,'\t%6.0f pts. trimmed, %6.0f s total',sum(sum(sum(~isnan(b)))),etime(clock,t1)),end
    end


    % -------------------------------------------------------------------
    % * write final images -- number of trimmed outliers image
    % -------------------------------------------------------------------
    if dofilt
        V.fname = ['ntrimmed.img'];
        V.descrip = ['Number of outliers trimmed from timeseries at each voxel.'];
        spm_write_vol(V,ntrimmed);

    end


    if vb > 0, fprintf(1,'%6.0f s Total, to dir: %s ',etime(clock,t1),pwd),end

    Pw = [];

    return






    % -------------------------------------------------------------------
    %
    %
    %
    %
    % * SUB-FUNCTIONS
    %
    %
    %
    %
    % -------------------------------------------------------------------


function [fbetas,ntrimmed] = process_slice(slicei,P,S,PSDXS,TR,V,Vall,F,varargin)
    %[a,b] = process_slice(slicei,P,S,PSDXS,TR,V,mask(:,:,slicei),smoothlen);

    global dographics
    global dosmooth
    global doysmooth
    global dohtw
    global dofilt
    global writefilt
    global dowritebetas
    global f1
    global f2
    global f3

    global whpredictors
    global contrasts
    global connames

    numframes = [1 1 1 1 1 1];      % should not be used in this version!! will return error if dosmooth is on.

    mask = []; if length(varargin) > 0, mask = varargin{1};, end
    smoothlen = 6; if length(varargin) > 1, smoothlen = varargin{2};, end

    global dims
    global cols

    % -------------------------------------------------------------------
    % * load the slice
    % -------------------------------------------------------------------
    fprintf(1,'\tLoading data >')
    et = clock;
    if ~isempty(mask) & ~(sum(sum(sum(mask))) > 0)
        % skip it
        fbetas = NaN * zeros([dims(1:2) cols]);
        ntrimmed = NaN;
        fprintf(1,'...Empty slice...')
        return
    else
        sl = timeseries_extract_slice(Vall,slicei);

        % load it, check values for middle slice only
        %if slicei == round(dims(3)./2),
        %    [ts,sl] = timeseries4(round(dims./2),P);
        %else,
        %    [ts,sl] = timeseries4(round(dims./2),P,[],1);
        %end
    end

    fprintf(1,'loaded in %3.2f s.',etime(clock,et))
    if ~isempty(mask), sl(:,:,1) = sl(:,:,1) .* mask;, end

    if dographics,
        if ishandle(f1), figure(f1);,else, f1 = figure; set(gcf,'Color','w');, end
        subplot 221; imagesc(mask);,colormap(gray); title('Mask'),colorbar('horiz'),drawnow;
        subplot 222; imagesc(sl(:,:,2));,colormap(gray); title('Slice timepoint 2'),colorbar('horiz'),drawnow;
        subplot 223; cla
        subplot 224; cla
    end


    % -------------------------------------------------------------------
    % * filter data and store in fsl
    % -------------------------------------------------------------------
    nsize = size(sl);
    fsl = zeros(nsize);     % filtered output slices -- initialize
    nsize = nsize(1:2);

    ntrimmed = NaN * zeros(nsize);
    wvox = find(sl(:,:,1) ~= 0 & ~isnan(sl(:,:,1)));
    [i,j] = ind2sub(size(sl(:,:,1)),wvox);

    if dofilt
        fprintf(1,'\n\tFiltering data for %3.0f voxels > ',length(i))


        et = clock;
        for k = 1:length(i)

            % first image usually off -- replace w/mean of first 10 images
            y = squeeze(sl(i(k),j(k),:)); y(1) = mean(y(min(length(y),10)));

            % filter
            y = S * y;

            % trim
            [y,ntrimmed(i(k),j(k))] = trimts(y,3,[]);
            fsl(i(k),j(k),:) = y;   % fsl is filtered slices

            if k == 1000, fprintf(1,'%3.0f s per 1000 vox.',etime(clock,et)), end
        end

        if dographics,
            if ishandle(f2), figure(f2);subplot(2,1,1);,else, f2 = tor_fig(2,1);, end, cla
            plot(scale(squeeze(sl(round(dims(1)./2),round(dims(2)./2),:))),'k','LineWidth',2)
            plot(scale(squeeze(fsl(round(dims(1)./2),round(dims(2)./2),:))),'r','LineWidth',2),
            legend({'Raw' 'Filtered'}); title('Standardized data from center voxel')
            subplot(2,1,2); cla;
        end

        clear sl

        if dographics,
            figure(f1); subplot 223; imagesc(ntrimmed);,colormap(gray); title('ntrimmed'),colorbar('horiz'),drawnow;
        end

    else
        fsl = sl;   % skip filtering
    end % if dofilt

    % -------------------------------------------------------------------
    % * write filtered images
    % -------------------------------------------------------------------

    emptyimg = zeros(dims);     % in case we need to create a new volume, used later as well

    if writefilt
        warning off
        fprintf(1,'\n\tWriting f (filtered) plane > ')
        et = clock;

        for img = 1:size(P,1)

            V(img) = V(1);
            [d,f,e]=fileparts(deblank(P(img,:)));
            V(img).fname = fullfile(d,['f' f e]);

            if img == 1, fprintf(1,'%s .',V(img).fname), end

            % create volume, if necessary
            if ~(exist(V(img).fname) == 2), spm_write_vol(V(img),emptyimg);,end

            spm_write_plane(V(img),fsl(:,:,img),slicei);

        end
        fprintf(1,'Done in %3.0f s\n',etime(clock,et))
        warning on
    end


    if ~isempty(PSDXS),

        % -------------------------------------------------------------------
        % * smooth timeseries data before fitting
        % -------------------------------------------------------------------
        if doysmooth

            fprintf(1,'\tSmoothing timeseries for %3.0f voxels > ',length(i))

            et = clock;
            for k = 1:length(i)

                fsl(i(k),j(k),:) = F * squeeze(fsl(i(k),j(k),:));
                if k == 1000, fprintf(1,'%3.0f s per 1000 vox.',etime(clock,et)), end

            end

            if dographics,
                if ishandle(f2), figure(f2);subplot(2,1,1);,else, f2 = tor_fig(2,1);, end
                plot(scale(squeeze(fsl(round(dims(1)./2),round(dims(2)./2),:))),'b','LineWidth',2),
                legend({'Raw' 'Filtered' 'Smoothed'}); title('Standardized data from center voxel')
            end

        end  % if doysmooth

        % -------------------------------------------------------------------
        % * fit the model to each nonzero voxel
        % -------------------------------------------------------------------
        betas = NaN * zeros([dims(1:2) cols]);
        cons = NaN * zeros([dims(1:2) size(contrasts,1)]);

        fprintf(1,'\n\tFitting model matrix and estimating contrasts > ')

        et = clock;
        for k = 1:length(i)
            % b = pinv(S * X) * (S * y) do PSDX * S first
            % equivalent to: betas(i(k),j(k),:) = PSDX * (S * squeeze(sl(i(k),j(k),:)));

            betas(i(k),j(k),:) = PSDXS * squeeze(fsl(i(k),j(k),:));   %for NaN replacement: PSDXS(:,~isnan(y)) * y(~isnan(y));

            if any(isnan(betas(i(k),j(k),:))), disp('Some betas are NaN for non-empty voxels.  This should not happen.'); keyboard, end

            % contrasts
            if ~isempty(contrasts)
                cons(i(k),j(k),:) = contrasts * squeeze(betas(i(k),j(k),whpredictors));
            end

            if k == 1000, fprintf(1,'%3.0f s per 1000 vox.',etime(clock,et)), end
        end

        if dographics,
            figure(f1); subplot 224; imagesc(betas(:,:,end));,colormap(gray); title('Last beta (intercept)'),colorbar('horiz'),drawnow;
            try, saveas(gcf,['slices_and_trimming_' num2str(slicei)],'tif'),catch, disp('Error saving tiff.'), end
            figure(f2); subplot(2,1,2); plot(squeeze(betas(round(dims(1)./2),round(dims(2)./2),:)),'k','LineWidth',2)
        end

    elseif (dosmooth | dohtw | dowritebetas)   %  & isempty(PSDXS)
        % If we have no model, but we want to do stuff with betas
        % Assume that we've entered Beta images in timeseries as slice
        betas = fsl;
    else   %  & isempty(PSDXS)
        % -------------------------------------------------------------------
        % * skip timeseries smooth and model fit if no model matrix is entered
        % -------------------------------------------------------------------
        % We're done, nothing left to do.
        betas = NaN * zeros([dims(1:2) cols]);
        return,
    end






    % -------------------------------------------------------------------
    % * smooth beta series within condition to regularize
    % -------------------------------------------------------------------
    st = cumsum([1 numframes]);
    en = st(2:end) - 1;         % ending values
    st = st(1:end-1);           % starting values

    flen = max(1,round(smoothlen/TR));    % smoothing filter length: 6 s to zero weight in exponential LP smooth

    if dosmooth
        fprintf(1,'\n\tRegularizing betas > ')
        et = clock;
        fbetas = NaN*zeros(size(betas));    %initialize to NaN to avoid dim problems

        for bet = 1:length(st)      % for each beta series

            % get LP filter
            [x,F] = smooth_timeseries(squeeze(betas(round(dims(1)./2),round(dims(1)./2),st(bet):en(bet))),flen);

            % loop through voxels, save in fbetas
            for k = 1:length(i)

                fbetas(i(k),j(k),st(bet):en(bet)) = F * squeeze(betas(i(k),j(k),st(bet):en(bet)));   %for NaN replacement: PSDXS(:,~isnan(y)) * y(~isnan(y));
                if k == 1000 & bet == 1, fprintf(1,'%3.0f s per 1000 vox.',etime(clock,et)), end

                if any(isnan(fbetas(i(k),j(k),st(bet):en(bet)))), disp('Smoothing has introduced NaN beta values.  This should not happen.'); keyboard, end

            end
        end

        if dographics,
            if ishandle(f2), figure(f2);subplot(2,1,1);,else, f2 = tor_fig(2,1);, end
            figure(f2); subplot(2,1,2); plot(squeeze(fbetas(round(dims(1)./2),round(dims(2)./2),:)),'r','LineWidth',2)
            legend({'Beta series' 'Smoothed'})
            saveas(gcf,['timeseries_center_voxel_' num2str(slicei)],'tif')
        end

        fprintf(1,'Done in %3.0f s\n',etime(clock,et))

    else
        % no smoothing
        fbetas = betas;

    end % if dosmooth

    clear betas


    % -------------------------------------------------------------------
    % * write betas (smoothing opt should be off!) only for preds to save
    % -------------------------------------------------------------------
    if dowritebetas
        fprintf(1,'\n\tWriting beta plane > ')
        et = clock;
        emptyimg = zeros(dims);     % in case we need to create a new volume

        bname = write_beta_slice(slicei,V(1),fbetas,emptyimg,'pmbeta_',whpredictors);

        fprintf(1,'%s. Done in %3.0f s\n',bname,etime(clock,et))
    end


    % -------------------------------------------------------------------
    % * write contrast values
    % -------------------------------------------------------------------
    if ~isempty(contrasts)

        if isempty(connames)
            for cc = 1:size(cons,3), connames{cc} = ['Contrast' num2str(cc)];,end
        end

        fprintf(1,'\n\tWriting contrast plane > ')
        et = clock;
        emptyimg = zeros(dims);     % in case we need to create a new volume

        cname = write_con_slice(slicei,V(1),cons,emptyimg,connames,contrasts);

        fprintf(1,'%s and other contrast imgs: Done in %3.0f s\n',cname{1},etime(clock,et))
    end


    % -------------------------------------------------------------------
    % * get height, width, delay parameters
    % -------------------------------------------------------------------
    if dohtw

        hconst = round(30 ./ TR);   % constrain peak height to within first 30 s
        hconst = min([repmat(hconst,1,length(numframes)); (numframes)]);  %min([hconst min(numframes)]);  % min or we get an error
        fprintf(1,'\n\tExtracting height, delay, width using first %3.0f betas > ',hconst)

        sz = size(fbetas); sz = [sz(1:2) length(st)];
        h = NaN*zeros(sz);    %initialize to NaN to avoid dim problems
        t = NaN*zeros(sz);    %initialize to NaN to avoid dim problems
        w = NaN*zeros(sz);    %initialize to NaN to avoid dim problems

        warning off
        et = clock;
        for bet = 1:length(st)      % for each beta series

            % loop through voxels, save in fbetas
            for k = 1:length(i)

                % use version 2, which parameterizes from first beta value
                [h(i(k),j(k),bet),t(i(k),j(k),bet),w(i(k),j(k),bet)] = fir2htw2(squeeze(fbetas(i(k),j(k),st(bet):en(bet))),hconst(bet),0,{'k'});

                if k == 1000 & bet == 1, fprintf(1,'%3.0f s per 1000 vox.',etime(clock,et)), end

            end
        end
        warning on

        if dographics,
            if ishandle(f3), figure(f3);,
            else, f3 = figure; set(gcf,'Color','w');,
                tmp = get(gcf,'Position'); tmp(3) = tmp(3)+.5*tmp(3);
                set(gcf,'Position',tmp)
            end
            subplot(1,4,1);
            imagesc(h(:,:,1));,colormap(hot); title('Height','FontSize',24),colorbar('horiz'),axis image; drawnow;
            subplot(1,4,2); imagesc(t(:,:,1));,title('Time to peak','FontSize',24),colorbar('horiz'),axis image,drawnow;
            subplot(1,4,3); imagesc(w(:,:,1));,title('Width','FontSize',24),colorbar('horiz'),axis image,drawnow;
            drawnow

            k = find(h == max(h(:)));    % max height
            [ii,jj,kk] = ind2sub(size(h),k);

            subplot(1,4,4),cla
            y = squeeze(fbetas(ii,jj,st(kk):en(kk)));
            plot(y,'ro-','LineWidth',2)
            [h2,t2,w2] = fir2htw2(y,hconst(kk),1,{'k'});

            saveas(gcf,['height_time_wid_' num2str(slicei)],'tif')
        end

        % write images for all beta series
        write_beta_slice(slicei,V(1),h,emptyimg,['height_event_']);

        write_beta_slice(slicei,V(1),t,emptyimg,['delay_event_']);

        write_beta_slice(slicei,V(1),w,emptyimg,['width_event_']);


        fprintf(1,'\tDone in %3.0f s\n',etime(clock,et))

    end     % do HTW


    return











function Pw = write_beta(bi,V,betas)
    % Volumetric version
    if bi < 10, myz = '000';, elseif bi < 100, myz = '00';, else myz = '000';,end
    V.fname = ['dx_beta_' myz num2str(bi) '.img'];
    V.descrip = ['Beta image for column ' num2str(bi)];

    spm_write_vol(V,betas);
    Pw = which(V.fname);

    return



function Pw = write_beta_slice(slicei,V,betas,emptyimg,varargin)
    % Pw = write_beta_slice(sliceindex,V,data,empty,prefix,whpredictors number list for naming)
    % Slice-o-matic version
    warning off % due to NaN to int16 zero conversions
    V.dim(4) = 16; % set to float to preserve decimals

    prefix = 'dx_beta_';
    if length(varargin) > 0, prefix = varargin{1};,end

    whpredictors = 1:size(betas,3);
    if length(varargin) > 1, whpredictors = varargin{2};,end

    for voli = whpredictors  % for each image/beta series point
        if voli < 10, myz = '000';, elseif voli < 100, myz = '00';, else myz = '000';,end
        V.fname = [prefix myz num2str(voli) '.img'];
        V.descrip = ['Beta image for column ' num2str(voli)];

        % create volume, if necessary
        if ~(exist(V.fname) == 2), spm_write_vol(V,emptyimg);,end

        spm_write_plane(V,betas(:,:,voli),slicei);
    end

    Pw = which(V.fname);

    warning on
    return



function Pw = write_con_slice(slicei,V,betas,emptyimg,varargin)
    % Pw = write_con_slice(sliceindex,V,data,empty,prefix,contrasts)
    % Slice-o-matic version
    % prefix is a cell array of names

    warning off % due to NaN to int16 zero conversions
    V.dim(4) = 16; % set to float to preserve decimals

    % default names
    for coni = 1:size(betas,3)
        prefix{coni} = ['Contrast_' num2str(coni)];
    end
    contrasts = zeros(size(betas,3),1);

    if length(varargin) > 0, prefix = varargin{1};,end
    if length(varargin) > 1, contrasts = varargin{2};,end

    for voli = 1:size(betas,3)  % for each image/beta series point

        V.fname = [prefix{voli} '.img'];
        V.descrip = ['Contrast image: ' num2str(contrasts(voli,:))];

        % create volume, if necessary
        if ~(exist(V.fname) == 2), spm_write_vol(V,emptyimg);,end

        spm_write_plane(V,betas(:,:,voli),slicei);

        Pw{voli} = which(V.fname);
    end



    warning on
    return
