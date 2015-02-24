function varargout=scn_roi_gui(Action,varargin)
%
%
% Tor Wager
% 
% Thanks to Tom Nichols for the excellent GUI shell!

%-----------------------------functions-called------------------------
%
%-----------------------------functions-called------------------------

% global variables we need for this shell

global cl
global EXPT
    
%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action='Init'; end


switch lower(Action)
    
    case lower('Init')
    %=======================================================================

    %clc
    %BrainVowager_defaults;
    scn_roi_gui('AsciiWelcome')
    scn_roi_gui('CreateMenuWin')

    % load EXPT

    if ~isempty(EXPT)
        disp('Using EXPT already in memory.'); 
    elseif exist('EXPT.mat')==2
        disp('loading EXPT.mat'); load EXPT; 
    else
        disp('You need to go to the directory containing individual subject results dirs ');
        disp('this directory should have information about the experiment stored in the file EXPT.mat')
        disp('No EXPT file in current directory.');
        EXPT = [];
        fprintf(1,'\n')
    end
        
    if ~isempty(cl)
        disp('Using clusters (cl variable) already in memory.'); 
    else
        disp('No clusters (cl variable) loaded yet.  May not be declared global in  base workspace?');
    end
    
    varargout{1} = EXPT;    
        
    case lower('AsciiWelcome')
    %=======================================================================
    disp( 'Welcome to scn_roi_gui: Extraction and display of ROI data.')
    fprintf('\n')

    case lower('Ver')
    %=======================================================================
    varargout = {'SCNlab Roi Menu'};
  


    case lower('CreateMenuWin')
    %=======================================================================
    close(findobj(get(0,'Children'),'Tag','scn_roi_gui Menu'))

    
    %-Initialize scn_roi_gui menu window
    %-----------------------------------------------------------------------
    [F, winwid, winh] = scn_roi_gui('initFigure');
    
    
    % default button sizes and positions, etc.
    
    topbutton = winh-100;        % y location of top button
    butspace = 35;               % y spacing of buttons
    
    fullbutxy = [165 30];       % full-length button width and height
    halfbutxy = [80 30];        % (left-hand) half-width button w and h
    rightbutx = 110+halfbutxy(1)+5;  % right-hand button start x
    
    



    %-Frames and text
    %-----------------------------------------------------------------------
    axes('Position',[0 0 80/winwid winh/winh],'Visible','Off')
    text(0.5,0.475,'Cluster and ROI tools',...
        'FontName','Times','FontSize',36,...
        'Rotation',90,...
        'VerticalAlignment','middle','HorizontalAlignment','center',...
        'Color',[1 1 1]*.6);

    text(0.2,0.96,'SCN Lab',...
        'FontName','Times','FontSize',16,'FontAngle','Italic',...
        'FontWeight','Bold',...
        'Color',[1 1 1]*.6);

    uicontrol(F,'Style','Frame','Position',[095 005 winwid-100 winh - 30],...
        'BackgroundColor',scn_roi_gui('Color'));  % colored frame
    uicontrol(F,'Style','Frame','Position',[105 015 winwid-120 winh - 50]);  % inner gray frame

    %-Buttons to launch scn_roi_gui functions
    %-----------------------------------------------------------------------
    uicontrol(F,'Style','Text',...
        'String','Get clusters','FontSize',14,...
        'HorizontalAlignment','Center',...
        'Position',[115 topbutton+30 fullbutxy],...
        'ForegroundColor','y','FontWeight','b');

    % -------------------------------------------
    % First section - Load and get
    % -------------------------------------------      
    
    % Create clusters from a mask
    spmg = 'P = spm_get(1,''*img'',''Select image.'');';
    str = [spmg 'cl = mask2clusters(P);'];                 % callback function
    str = scn_roi_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Create(mask)',...
        'Position',[110 topbutton-(1-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    

    % Load an existing clusters .mat file
    spmg = 'P = spm_get(1,''*mat'',''Select mat file containing cl variable (clusters).'');';
    str = [spmg ',load(P);'];                 % callback function
    str = scn_roi_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Load(.mat)',...
        'Position',[rightbutx topbutton-(1-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
   
    % Draw cluster ROI mask
    str1 = ['EXPT = CreateExpt(''CreateField'',''mask''); '];
    str1b = ['disp([''Masking spheres with image: '' EXPT.mask]);'];
    str1c = ['if ~isfield(EXPT,''overlay''), EXPT.overlay=[];, disp([''Overlay is: '' EXPT.overlay]);'];
    str2 = ['bilat = input(''Make ROI spheres bilateral?''); '];
    str3 = ['cl = sphere_roi_tool(''mask'',EXPT.mask,''bilat'',bilat,''overlay'',EXPT.overlay); '];
    
    str = [str1 str1b str2 str3];                 % callback function
    str = scn_roi_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Draw(spheres)',...
        'Position',[110 topbutton-(2-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    % get clusters from Talairach
    str = ['talairach_cluster_plugin;'];                  % callback function
    str = scn_roi_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Talairach',...
        'Position',[rightbutx topbutton-(2-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
        
    % -------------------------------------------
    % Next section - Select and parcel
    % -------------------------------------------   
    
    uicontrol(F,'Style','Text',...
    'String','Prune and parcellate','FontSize',14,...
    'HorizontalAlignment','Center',...
    'Position',[115 topbutton-(3-1)*butspace fullbutxy],...
    'ForegroundColor','y','FontWeight','b');
        
        
    % parcel clusters based on anatomy (hierarchical clustering).
    str = ['cl = anat_subclusters(cl);'];                  % callback function
    str = scn_roi_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Parcel(anat)',...
        'Position',[rightbutx topbutton-(4-1)*butspace+15 halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
 
    % parcel clusters based on functional data (k-means clustering).
    str1 = ['CLU = clusters2CLU(cl);'];
    str2 = ['[cl,nclasses,colors] = cluster_kmeans_parcel(x,CLU,1);'];
    
    str = [str1 str2];                                      % callback function
    str = scn_roi_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Parcel(func)',...
        'Position',[110 topbutton-(4-1)*butspace+15 halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');


    % Graphic select: save in cl only after you press 'x'!
    str1 = ['CLU = clusters2CLU(cl);'];
    str2 = ['[clout,cl] = cluster_graphic_select(cl);'];
    
    str = [str1 str2];                                      % callback function
    str = scn_roi_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','GraphicSelect',...
        'Position',[110 topbutton-(5-1)*butspace+15 fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b')
    
        
    % -------------------------------------------
    % Next section - View
    % -------------------------------------------   
    
%     uicontrol(F,'Style','Text',...
%     'String','View','FontSize',14,...
%     'HorizontalAlignment','Center',...
%     'Position',[115 topbutton-(6-1)*butspace fullbutxy],...
%     'ForegroundColor','y','FontWeight','b');

    % -------------------------------------------

    % Orthviews Menu
    buttontext = {'Viewer Menu' 'New Orthviews' 'Add to orthviews' 'Orth in unique colors' 'Axial Montage' 'Medial Montage' 'Red/Blue Orthviews'};
        
    str = 'scn_roi_gui(''viewer'');';                          % callback function
    pop1 = uicontrol(F,'Style','popupmenu', 'Tag','ViewerPop', ...
    'String',buttontext,...
        'Position',[110 topbutton-(6-1)*butspace+15 fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
  
    buttontext = {'red' 'green' 'blue' 'yellow' 'orange' 'cyan' 'magenta'};    
    str = [];                          % callback function
    pop1 = uicontrol(F,'Style','popupmenu', 'Tag','ColorPop', ...
    'String',buttontext,...
        'Position',[110 topbutton-(7-1)*butspace+15 fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    buttontext = {'Overlay Image' 'Single subj' 'EXPT.overlay' 'Choose your own'};
    %which('scalped_single_subj_T1.img');    
    str = 'scn_roi_gui(''overlay'');';                          % callback function
    pop1 = uicontrol(F,'Style','popupmenu', 'Tag','OverlayPop', ...
    'String',buttontext,...
        'Position',[110 topbutton-(8-1)*butspace+15 fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');      
        

 
    % -------------------------------------------
    % Get data extraction menu
    str = 'scn_extract_data_gui;';
    str = scn_roi_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Extract Data Menu','FontSize',14,...
        'Position',[110 topbutton-(9-1)*butspace fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','y','FontWeight','b');
        
    % -------------------------------------------
    
  
    



    set(F,'Pointer','Arrow','Visible','on')



    case lower('Color')
    %=======================================================================
    % scn_roi_gui('Color')
    %-----------------------------------------------------------------------
    % %-Developmental livery
    % varargout = {[0.7,1.0,0.7], 'Lime Green'};
    %-Distribution livery
    varargout = {[.5 0.7 .5], 'Purple'};

    
    case lower('ExpandString')
    %=======================================================================
    % scn_roi_gui('ExpandString')
    % Expand an action button callback string (a command string to be
    % evaluated)
    % so that it first displays the command, and then executes it
    %-----------------------------------------------------------------------      
        str = varargin{1}; str2 = [];
        for i = 1:length(str)
            if str(i) == '''', str2(end+1) = ''''; str2(end+1) = '''';
            else, str2(end+1) = str(i);
            end
        end

        str = ['disp(''' char(str2) '''), ' str ];     % display then execute

        varargout = {str};
  
        
    case lower('initFigure')
    %=======================================================================
    % [F, winwid, winh] = scn_roi_gui('initFigure')
    %-----------------------------------------------------------------------
    % Get the position of the main BrainVowager menu, or if
    % not available, default screen pos.
        
        % default sizes, etc.
    S = get(0,'ScreenSize');

    winwid = 300;               % window width
    winh = 400;                 % window height
    pos = [S(3)/2+150,S(4)/2-140,winwid,winh];  % default
    
    % Look for existing figures in order:
    h = [];

    h = findobj('Tag','robust_toolbox_gui Menu');
    h = findobj('Tag','scn_fir_results_gui Menu');
    h = findobj('Tag','BrainVowager_gui Menu');,
    h = findobj('Tag','Meta_Analysis_gui Menu');,
    
    if ~isempty(h), 
        pos = get(h,'Position');, 
        winwid = pos(3); winh = pos(4);
        pos(2) = pos(2) - winh;   % put next to main figure
    end
    
    %-Open scn_roi_gui menu window
    %----------------------------------------------------------------------
    
    F = figure('Color',[1 1 1]*.8,...
        'Name',scn_roi_gui('Ver'),...
        'NumberTitle','off',...
        'Position',pos,...
        'Resize','off',...
        'Tag','scn_roi_gui Menu',...
        'Pointer','Watch',...
        'MenuBar','none',...
        'Visible','off');
    
    varargout{1} = F; varargout{2} = winwid; varargout{3} = winh;
    

    case lower('viewer')
    %=======================================================================
    % scn_roi_gui('viewer')
    %----------------------------------------------------------------------
    
    global cl
    if isempty(cl), warning('cl is empty.  not created, or not declared global in base workspace?');, end
    
    % get overlay image
    scn_roi_gui('overlay');
    han = findobj('Tag','scn_roi_gui Menu');
    P = guidata(han);
    
    % callbacks
    %buttontext = {'viewer menu' 'New Orthviews' 'Add to orthviews' 'Orth in unique colors' 'Axial Montage' 'Medial Montage'};
    pop1 = findobj('Tag','ViewerPop'); indx = get(pop1,'Value');
    
    % colors
    %buttontext = {'red' 'green' 'blue' 'yellow' 'orange' 'cyan' 'magenta'};    
    pop1 = findobj('Tag','ColorPop'); cindx = get(pop1,'Value');
    

    colors1 = {'r' 'g' 'b' 'y' 'r' 'c' 'm'};
    colors2 = {[1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 .5 0] [0 1 1] [1 0 1]};  
    
    % for montages
    if indx == 4 | indx == 5, 
        color = colors1{cindx};
    else
        % for others
        color = colors2{cindx};
    end

    % callback functions
    
    callbk = [{'0;'} ...
            {'cluster_orthviews(cl,{color},''overlay'',P);'} ...
           {'cluster_orthviews(cl,{color},''add'');'} ...
           {'cluster_orthviews(cl,''unique'',''overlay'',P);'} ...
           {'montage_clusters(P,cl,{color});'} ...
           {'montage_clusters_medial(P,cl,{color});'} ...
           {'cluster_orthviews(cl,''bivalent'',''overlay'',P);'} ...
           ];
       
    eval(callbk{indx});
    
    case lower('overlay')
    %=======================================================================
    % scn_roi_gui('overlay')
    % choose overlay image for montages, etc.
    %----------------------------------------------------------------------
    
    pop1 = findobj('Tag','OverlayPop'); 
    
    if isempty(pop1), disp('For overlay selection, try scn_roi_gui'); indx = 3;, 
    
    else
        indx = get(pop1,'Value');
    end
    
    %buttontext = {'Overlay Image' 'Single subj' 'EXPT.overlay' 'Choose your own'};
    
    switch indx
        case 1, P = which('scalped_single_subj_T1.img'); 
        case 2, P = which('scalped_single_subj_T1.img'); 
        case 3, P = which(EXPT.overlay); 
        case 4, P = spm_get(1,'*IMAGE','Choose overlay image.');
    end
    
    if isempty(P), warning('Cannot find overlay image on path!');, end
    
    % add to gui data
    han = findobj('Tag','scn_roi_gui Menu');
    guidata(han,P);
    
    otherwise
    %=======================================================================
    error('Unknown action string')

    %======================================================================
    
end


return


