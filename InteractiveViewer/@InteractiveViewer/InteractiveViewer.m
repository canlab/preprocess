% InteractiveViewer is a tool for displaying a reference brain with SPM and monitoring
% the cursor position to display auxiliary information, usually from an associated data set.
% Examples include displaying a voxel's timecourse over a session, or displaying a list of
% nearby clusters from a meta-analysis.
%
% InteractiveViewer() creates a default InteractiveViewer and initializes the SPM Graphics window, if absent
%
% InteractiveViewer(['DisplayVol', display_img | 'AssociatedVols', data_imgs | 'IVObserver', observer | 'UseExistingGraphicsWindow', [0|1] | 'LoadDataOnDemand', [0|1] ])
%   Initializes the SPM window with the display image, sets up the data, and attaches the observer. If you enter 
%   'UseExistingGraphicsWindow', 1, it will *not* initialize a new window. Multiple observers may be specified. 
%   Alternatively, the observers can be attached later. NB: They are called "observers" because in the Observer 
%   pattern, they "watch" for changes somewhere else, in this case, the currently selected voxel location.
%
%   display_img: a filename or an spm_vol structure
%   data_imgs: a char matrix of filenames, a cellstr of filenames, an array of spm_vol structure, or the most useful form,
%           a structure with two fields, 'name' and 'data_imgs' - 'name' contains the name of the group, and 'data_imgs' contains
%           one of three earlier possibilities for data_imgs
%   observer: an object of class IVObserver. There are several existing functions for generating common observers, such as 
%           errorBarPlotObserver(), tsPlotObserver(), and more. See the InteractiveViewer directory for more examples.
%
% The IVObservers receive a reference to the InteractiveViewer, which they should query to retrieve data
% associated with the InteractiveViewer. Common forms of data are queryable, but something truly unique to
% a particular IVObserver should not be placed in the InteractiveViewer. All IVObservers receive the mm and
% voxel position of the cursor for the display image (the thing they *observe*). 

% Here are the queryable properties of IV:
%   'CurrentVols': the spm_vol structures of the associated data imgs
%	'CurrentVolsData': the 4-D values of the associated data imgs
%   'CurrentVolsVoxPos': the voxel position in the space of the data imgs, not the display img
%	'CurrentVolsFilenames': the filenames of the data imgs
%	'CurrentVolsTs': the timeseries of the voxel data across all data imgs
%   'NumDataSets', 'NumDataVols': the number of image sets
%   'CurrentDataIdx', 'CurrentDataVolsIdx': Which data set is the current one
%   'DisplayVol': spm_vol structure of the image being overlaid upon
%   'LoadDataOnDemand': boolean as to whether or not to load data for images all at once or on demand (default: 0)
%   'IdenticalSpace': boolean as to whether or not the data sets should be in the same image space (default: 0)
%   Data set properties:
%       'Vols': spm_vol structures of the data set
%       'Data': all the data loaded so far of the set (WARNING: this could be a LOT of data if LoadDataOnDemand is off)
%       'Ts', 'Timeseries': a vector of the data across all images in the set at the current cursor position
%       'VoxPos', 'CurrentVoxPos': current voxel position for the data set
%       'Filenames': list of filenames of the data set
%       
%
% InteractiveViewer has no copy constructor, array or cell referencing, as it is a Singleton. It
% has been constructed as a reference object, hence, all copies share the same data. This
% is due to SPM keeping only one Graphics window (as of SPM2), and avoids issues with conflicting
% callbacks. InteractiveViewer is also the Subject in an Observer pattern. For more on the Observer
% and Singleton patterns, see:
%   Erich Gamma, Richard Helm, Ralph Johnson, and John Vlissides.
%   Design Patterns: Elements of Reusable Object-Oriented Software. Addison-Wesley, 1995.
%
% Created by Matthew Davidson, 2006
%
% Examples:
% to create a new iv (interactive viewer) object:
% iv = InteractiveViewer('UseExistingGraphicsWindow', 1, 'AssociatedVols', images);
%
% to attach a data viewing pane to the iv object:
% obs3 = tsAncovaPlotObserver(plotAgainst, controlFor, groupBy);
% attach(iv, obs3);
%
% to get data from an existing iv:
% ts = getDataVolProp(iv, 1, 'Ts');
%
% get data from an iv with multi-subject data:
% vox_data = cell(1, length(X)); % X is cell array of X data in mediation
% analysis, for example.
% for i = 1:length(X), vox_data{i} = getDataVolProp(iv, i, 'Ts'); end

function iv = InteractiveViewer(varargin)
    useExistingGraphicsWindow = 0;
    loadDataOnDemand = 0;
    identicalSpace = 0;
    
    hSpmFig = [];
    hPlotFig = [];
    hControlPanel = [];
    hPlotPanel = [];
    hPlotAxes = [];
    hCurrentDataVolsPopup = [];
    

    currentMmPos = [0 0 0]';
    displayImg = struct('vol', [], 'currentVoxPos', []);

    dataVols = struct('vol', {}, 'currentVoxPos', {}, 'name', {}, 'data', {});
    currentDataVolsIdx = 0;

    clusters = struct([]);
    currentClusterIdx = 0;
    
    observers = struct('IVObs', {}, 'name', {}, 'attached', {});
    observers_to_attach = IVObserver([]);

    predictors = struct('data', {}, 'name', {});

    iv.attach = @attach_;
    iv.detach = @detach_;
    iv.notify = @notify_;
    iv.display = @display_;
    iv.get = @get_;
    iv.getDataVolProp = @getDataVolProp_;
    iv.set = @set_;
    iv = class(iv, 'InteractiveViewer');

    parseInputs(varargin);
    setupSpmWindow();
    setupPlotWindow();

    computeDisplayImgVoxPos();
    if(~isempty(dataVols))
        computeDataVoxPos();
    end

%     for i=1:length(observers_to_attach)
%         attach_(observers_to_attach(i));
%     end
    attach_(observers_to_attach);


    %%%%%%%%%%%%%%%%%%%%%%
    % Accessor functions
    %%%%%%%%%%%%%%%%%%%%%%

    function val = get_(propName)
        switch propName
            case {'NumDataSets', 'NumDataVols'}
                val = length(dataVols);
            case 'CurrentVols'
                val = getDataVolProp_(currentDataVolsIdx, 'Vols');
            case 'CurrentVolsData'
                val = getDataVolProp_(currentDataVolsIdx, 'Data');
            case 'CurrentVolsVoxPos'
                val = getDataVolProp_(currentDataVolsIdx, 'VoxPos');
            case 'CurrentVolsFilenames'
                val = getDataVolProp_(currentDataVolsIdx, 'Filenames');
            case 'CurrentVolsTs'
                val = getDataVolProp_(currentDataVolsIdx, 'Ts');
            case 'DisplayVol'
                val = displayImg.vol;
            case 'NumObservers'
                val = length(observers);
            case 'UseExistingGraphicsWindow'
                val = useExistingGraphicsWindow;
            case 'LoadDataOnDemand'
                val = loadDataOnDemand;
            case {'CurrentDataIdx', 'CurrentDataVolsIdx'}
                val = currentDataVolsIdx;
            otherwise
                error([propName ' is not a valid property']);
        end
    end

    function val = getDataVolProp_(idx, propName)
        if(~isscalar(idx))
            error('idx must be a single scalar in %s\n', mfilename());
        else
            switch propName
                case 'Vols'
                    val = dataVols(idx).vol;
                case 'Data'
                    if(loadDataOnDemand)
                        val = spm_read_vols(dataVol.vol);
                    else
                        val = dataVols(idx).data;
                    end
                case {'Ts', 'Timeseries'}
                    dataVoxPos = dataVols(idx).currentVoxPos;
                    if(loadDataOnDemand)
                        val = arrayfun(@(currentVol) spm_sample_vol(currentVol, dataVoxPos(1), dataVoxPos(2), dataVoxPos(3), 0), dataVols(idx).vol);
                    else
                        val = squeeze(dataVols(idx).data(dataVoxPos(1), dataVoxPos(2), dataVoxPos(3), :));
                    end
                case {'VoxPos', 'CurrentVoxPos'}
                    val = dataVols(idx).currentVoxPos;
                case 'Filenames'
                    val = {dataVols(idx).vol.fname};
                otherwise
                    error([propName ' is not a valid property']);
            end
        end
    end

    function set_(varargin)
        for i=1:2:length(varargin)
            propName = varargin{i};
            val = varargin{i+1};
            switch propName
                case 'CurrentVols'
                    dataVols(currentDataVolsIdx) = makeDataVol(spmVol(val));
                    computeDataVoxPos(currentDataVolsIdx);
                case 'DisplayVol'
                    displayImg.vol = spmVol(val);
                    spm_image('init', displayImg.vol);
                case {'CurrentDataIdx', 'CurrentDataVolsIdx'}
                    if(val > length(dataVols) || val < 1)
                        error('There are only %d img data sets. Cannot set primary data set to %d\n', length(dataVols), val);
                    elseif(currentDataVolsIdx ~= val)
                        currentDataVolsIdx = val;
                        notify_();
                    end
                otherwise
                    error([propName ' is not a valid property']);
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Observer pattern functions
    %%%%%%%%%%%%%%%%%%%%%%
    
    function attach_(obs)
        if(~isa(obs, 'IVObserver'))
            error('obs is not of class (or a subclass of) IVObserver');
        else
            oldNumObservers = length(observers);
            for i=1:length(obs)
                obs(i) = setIv(obs(i), iv);
                observers(end+1) = struct('IVObs', obs(i), 'name', [], 'attached', 1);
            end
            notify_(oldNumObservers+1:oldNumObservers+length(obs));
        end
    end

    function detach_(idxs)
        observers(idxs) = [];
        notify_();
    end

    function notify_(idxs)
        if(~isempty(dataVols))
            if(~exist('idxs', 'var') || isempty(idxs))
                idxs = 1:length(observers);
            end

            computeDataVoxPos();

            set(0, 'CurrentFigure', hPlotFig);
            axes('Parent', hPlotPanel, 'Visible', 'off');
            for i=idxs
                subplot(length(observers), 1, i, 'align');
                update(observers(i).IVObs, currentMmPos, displayImg.currentVoxPos);
            end
        end
    end

    function display_()
        display(observers);
    end
    

    %%%%%%%%%%%%%%%%%%%%%%
    % Init functions
    %%%%%%%%%%%%%%%%%%%%%%
    
    function parseInputs(args)
        DEFAULT_DISPLAY_IMG = 'single_subj_T1.img';

        if(~isempty(args) && isa(args{1}, 'InteractiveViewer'))
            error('Only one InteractiveViewer object is allowed.');
        end
        
        %% First pass - needed before handling certain other params
        for i=1:2:length(args)
            if(ischar(args{i}))
                switch(args{i})
                    case 'DisplayVol'
                        displayImg.vol = spmVol(args{i+1});
                    case 'IdenticalSpace'
                        identicalSpace = args{i+1};
                    case 'LoadDataOnDemand'
                        loadDataOnDemand = args{i+1};
                    case 'UseExistingGraphicsWindow'
                        useExistingGraphicsWindow = args{i+1};
                end
            end
        end

        %% Second pass
        for i=1:2:length(args)
            if(ischar(args{i}))
                switch(args{i})
                    case {'DataSet' 'AssociatedVols' 'AssociatedDataVols' 'AssociatedDataSet'}
                        if(isstruct(args{i+1}) && isfield(args{i+1}, 'name') && isfield(args{i+1}, 'data_imgs'))
                            newDataVol = makeDataVol(spmVol(args{i+1}.data_imgs));
                            newDataVol.name = args{i+1}.name;
                        else
                            newDataVol = makeDataVol(spmVol(args{i+1}));
                        end
                        addDataVols(newDataVol);
                    case 'Cluster'
                        addCluster(args{i+1});
                    case 'IVObserver'
                        observers_to_attach = append(observers_to_attach, args{i+1});
                        
                    case {'UseExistingGraphicsWindow', 'LoadDataOnDemand', 'IdenticalSpace', 'DisplayVol'}
                        % do nothing
                        
                    otherwise
                        warning('InteractiveViewer:unknownKeywordParameter', 'Unrecognized keyword parameter: "%s"', args{i});
                end
            end
        end
        
        if(identicalSpace)
            verifyIdenticalSpace(dataVols);
        end
        
        if(isempty(displayImg.vol))
            displayImg.vol = spmVol(which(DEFAULT_DISPLAY_IMG));
        end
    end
    
    function addDataVols(newDataVol)
        dataVols(end+1) = newDataVol;
        if(currentDataVolsIdx == 0)
            currentDataVolsIdx = 1;
        end
    end

    function addCluster(cl)
        if(isempty(clusters))
            clusters = cl;
        else
            clusters(end+1) = cl;
        end
        if(currentClusterIdx == 0)
            currentClusterIdx = 1;
        end
    end

    function setupSpmWindow()
        if(~useExistingGraphicsWindow)
            spm_image('init', displayImg.vol);
        end
        hSpmFig = spm_figure('GetWin', 'Graphics');

        if(~isempty(clusters))
            colors = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            for i=1:length(clusters)
                spm_orthviews('AddColouredBlobs', 1, clusters(i).XYZmm, clusters(i).Z, clusters(i).M, colors(i,:));
            end
            spm_orthviews('Reposition', clusters(1).center);
        end
        
        hReg = uicontrol(hSpmFig, 'Style', 'Text', 'String', 'InteractiveViewer hReg', ...
            'Position', [100 200 100 025], 'Visible', 'Off', ...
            'FontName', 'Times', 'FontSize', 14, 'FontWeight', 'Bold', ...
            'HorizontalAlignment', 'Center');
        hReg = spm_XYZreg('InitReg', hReg, displayImg.vol.mat, displayImg.vol.dim(1:3)');
        spm_XYZreg('Add2Reg', hReg, 0, @ivCallback);
        spm_orthviews('Register', hReg);
    end

    function ivCallback(command, mmpos, h, hReg) %#ok # mlint notice
        currentMmPos = mmpos;
        computeDisplayImgVoxPos();
        %fprintf(' mm pos: %0.3f %0.3f %0.3f\n', currentMmPos(1), currentMmPos(2), currentMmPos(3));
        %fprintf('Current vox pos: %0.3f %0.3f %0.3f\n', displayImg.currentVoxPos(1), displayImg.currentVoxPos(2), displayImg.currentVoxPos(3));
        notify_();
    end

    function setupPlotWindow()
        PLOT_TAG = 'InteractiveViewer Plot Window';
        BG_COLOR = [.8 .8 .8];
        close(findobj('Tag', PLOT_TAG));
        hPlotFig = figure('Name', 'InteractiveViewer', 'Tag', PLOT_TAG, 'Color', BG_COLOR);
        hControlPanel = uipanel('Parent', hPlotFig, 'Position', [0 .8 1 .2], 'BackgroundColor', BG_COLOR);
        hPlotPanel = uipanel('Parent', hPlotFig, 'Position', [0 0 1 .8], 'BackgroundColor', BG_COLOR);
        
        hCurrentDataVolsPopup = uicontrol('parent', hControlPanel, 'style', 'popupmenu', 'String', {dataVols.name}, 'Callback', @(hObject,eventdata)set_('CurrentDataIdx', get(hObject,'Value')));
    end
    
    function computeDisplayImgVoxPos()
        displayImg.currentVoxPos = round(mm2voxel(currentMmPos, displayImg.vol.mat));
    end

    function computeDataVoxPos(idxs)
        if(~isempty(dataVols))
            if(~exist('idxs', 'var') || isempty(idxs))
                idxs = 1:length(dataVols);
            end
            for i=idxs
                dataVols(i).currentVoxPos = round(mm2voxel(currentMmPos, dataVols(i).vol(1).mat));
            end
        end
    end
    
    function dataVol = makeDataVol(spmVol)
        dataVol = struct('vol', spmVol, 'currentVoxPos', [], 'name', spmVol(1).descrip, 'data', []);
        if(~loadDataOnDemand)
            dataVol.data = spm_read_vols(dataVol.vol);
        end
    end
end


%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%

function verifyIdenticalSpace(dataVols)
   allvols = vertcat(dataVols(:).vol);
   
   dims = vertcat(allvols.dim); 
   dims = dims(:,1:3);
   if(size(unique(dims, 'rows'), 1) > 1)
       error('Not all images have the same dimensions when they should.');
   end
   
   wh_mat_diff = find(squeeze(any(any(abs(diff(cat(3, allvols.mat), 1, 3)) > .0001))));
   if(~isempty(wh_mat_diff))
       for i=1:length(wh_mat_diff)
           fprintf('Image %s differs from %s\n', allvols(i).fname, allvols(i+1).fname);
       end
       error('Images required to be in identical images space, but aren''t.\n');
   end
end