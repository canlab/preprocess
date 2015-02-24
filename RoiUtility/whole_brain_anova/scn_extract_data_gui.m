function varargout=scn_extract_data_gui(Action,varargin)
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
    scn_extract_data_gui('AsciiWelcome')
    scn_extract_data_gui('CreateMenuWin')

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
    end
    
    varargout{1} = EXPT;    
        
    case lower('AsciiWelcome')
    %=======================================================================
    disp( 'Welcome to scn_extract_data_gui: Extraction and display of ROI data.')
    fprintf('\n')

    case lower('Ver')
    %=======================================================================
    varargout = {'SCNlab Roi Menu'};
  


    case lower('CreateMenuWin')
    %=======================================================================
    close(findobj(get(0,'Children'),'Tag','scn_extract_data_gui Menu'))

    
    %-Initialize scn_extract_data_gui menu window
    %-----------------------------------------------------------------------
    [F, winwid, winh] = scn_extract_data_gui('initFigure');
    
    
    % default button sizes and positions, etc.
    
    topbutton = winh-100;        % y location of top button
    butspace = 35;               % y spacing of buttons
    
    fullbutxy = [165 30];       % full-length button width and height
    halfbutxy = [80 30];        % (left-hand) half-width button w and h
    rightbutx = 110+halfbutxy(1)+5;  % right-hand button start x
    
    



    %-Frames and text
    %-----------------------------------------------------------------------
    axes('Position',[0 0 80/winwid winh/winh],'Visible','Off')
    text(0.5,0.475,'Data Extraction Menu',...
        'FontName','Times','FontSize',36,...
        'Rotation',90,...
        'VerticalAlignment','middle','HorizontalAlignment','center',...
        'Color',[1 1 1]*.6);

    text(0.2,0.96,'SCN Lab',...
        'FontName','Times','FontSize',16,'FontAngle','Italic',...
        'FontWeight','Bold',...
        'Color',[1 1 1]*.6);

    uicontrol(F,'Style','Frame','Position',[095 005 winwid-100 winh - 30],...
        'BackgroundColor',scn_extract_data_gui('Color'));  % colored frame
    uicontrol(F,'Style','Frame','Position',[105 015 winwid-120 winh - 50]);  % inner gray frame

    %-Buttons to launch scn_extract_data_gui functions
    %-----------------------------------------------------------------------
    uicontrol(F,'Style','Text',...
        'String','Contrast data','FontSize',14,...
        'HorizontalAlignment','Center',...
        'Position',[115 topbutton+30 fullbutxy],...
        'ForegroundColor','y','FontWeight','b');

    
    % -------------------------------------------
    % First section - Contrast data
    % -------------------------------------------      
    
    str = 'cl = extract_contrast_data(EXPT.SNPM.P,cl);';            % callback function    
    str = scn_extract_data_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Contrast',...
        'Position',[110 topbutton-(1-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    str = 'disp(''in progress.'';';                                 % callback function
    str = scn_extract_data_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Plot',...
        'Position',[rightbutx topbutton-(1-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');

        
  str = 'scn_extract_data_gui(''setup_seeds'');';            % callback function    
    str = scn_extract_data_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Get Seeds',...
        'Position',[110 topbutton-(2-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');        
        
        


    
    % -------------------------------------------
    % 2nd section - HRFs
    % -------------------------------------------      
    uicontrol(F,'Style','Text',...
        'String','HRF estimates','FontSize',14,...
        'HorizontalAlignment','Center',...
        'Position',[115 topbutton-(3-1)*butspace fullbutxy],...
        'ForegroundColor','y','FontWeight','b');
    
    
    
    str = '[cl,EXPT] = extract_dxbeta_data(EXPT,cl);';            % callback function
    str = scn_extract_data_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','FIR dx*',...
        'Position',[110 topbutton-(4-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');

    str = 'plot_dx_hrfs(EXPT,cl,0,1,EXPT.FIR.smoothlen,EXPT.FIR.indiv);';
    str = scn_extract_data_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','FIR dx*',...
        'Position',[rightbutx topbutton-(4-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');   
    
    
    %Batch time-course
    %str = ['EXPT = CreateExpt(''CreateField'',''dxnames'',''FIR'');', ...   
    %    'EXPT = robust_batch_dx_plot(EXPT,''both'');' ];       % callback function

    str = 'EXPT = robust_batch_dx_plot(EXPT,''both'');';
    str = scn_extract_data_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Batch brain FIR extract',...
        'Position',[110 topbutton-(5-1)*butspace fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    % -------------------------------------------

    
    % -------------------------------------------
    % Next section
    % -------------------------------------------
    
    uicontrol(F,'Style','Text',...
    'String','Raw data','FontSize',14,...
    'HorizontalAlignment','Center',...
    'Position',[115 topbutton-(6-1)*butspace fullbutxy],...
    'ForegroundColor','y','FontWeight','b');

    % -------------------------------------------
    
    str = 'cl = extract_raw_data(EXPT,cl);';
    str = scn_extract_data_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Extract raw data',...
        'Position',[110 topbutton-(7-1)*butspace fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
 
    set(F,'Pointer','Arrow','Visible','on')



    case lower('Color')
    %=======================================================================
    % scn_extract_data_gui('Color')
    %-----------------------------------------------------------------------
    % %-Developmental livery
    % varargout = {[0.7,1.0,0.7], 'Lime Green'};
    %-Distribution livery
    varargout = {[.5 0.7 .5], 'Purple'};

    
    case lower('ExpandString')
    %=======================================================================
    % scn_extract_data_gui('ExpandString')
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
    % [F, winwid, winh] = scn_extract_data_gui('initFigure')
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
    h = findobj('Tag','BrainVowager_gui Menu');,
    h = findobj('Tag','scn_roi_gui Menu');,
    
    if ~isempty(h), 
        pos = get(h,'Position');, 
        winwid = pos(3); winh = pos(4);
        pos(1) = pos(1) + winwid;   % put next to main figure
    end
    
    
    %-Open scn_extract_data_gui menu window
    %----------------------------------------------------------------------
    
    F = figure('Color',[1 1 1]*.8,...
        'Name',scn_extract_data_gui('Ver'),...
        'NumberTitle','off',...
        'Position',pos,...
        'Resize','off',...
        'Tag','scn_extract_data_gui Menu',...
        'Pointer','Watch',...
        'MenuBar','none',...
        'Visible','off');
    
    varargout{1} = F; varargout{2} = winwid; varargout{3} = winh;
    
    
        
    case lower('setup_seeds')
    %=======================================================================
    % cl = scn_extract_data_gui('setup_seeds')
    %-----------------------------------------------------------------------
    % Save data in cl.CONTRAST in EXPT.seeds for seed analyses with robseed
            
        EXPT = CreateExpt('robfit');
        fprintf(1,'\n');
        disp(EXPT.SNPM.connames)
        fprintf(1,'\n');
        wh = input('Which contrast? ');
        
        fprintf(1,'Using data from this contrast as seed data:\n%s\n',EXPT.SNPM.connames(wh,:));
        
        if isempty(cl), error('cl is empty!  not declared global in base workspace?');,end
        
        if ~isfield(cl,'CONTRAST')
            cl = extract_contrast_data(EXPT.SNPM.P,cl);
        end
        
        EXPT.seeds = [];
        for i = 1:length(cl)
            EXPT.seeds = [EXPT.seeds cl(i).CONTRAST.data(:,wh)];
        end
        
        % name seed clusters, if we don't have them already
        if ~isfield(cl,'shorttitle'),
            cl = cluster_names(cl);
        end
            
        for i= 1:length(cl), EXPT.seednames{i} = cl(i).shorttitle;, end
        disp('Seeds saved in EXPT.seeds');
        
    
    otherwise
    %=======================================================================
    error('Unknown action string')

    %======================================================================
    
end


return


