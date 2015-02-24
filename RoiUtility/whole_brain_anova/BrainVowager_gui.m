function varargout=BrainVowager_gui(Action,varargin)
%
%
% Tor Wager
% 
% Thanks to Tom Nichols for the excellent GUI shell!

%-----------------------------functions-called------------------------
%
%-----------------------------functions-called------------------------


%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action='Init'; end


switch lower(Action)
    
    case lower('Init')
    %======================================================================
    %=
    global EXPT
    clc
    %BrainVowager_gui_defaults;
    BrainVowager_gui('AsciiWelcome')
    BrainVowager_gui('CreateMenuWin')

 
    if ~isempty(EXPT)
        disp('Using EXPT already in memory.'); 
    elseif exist('EXPT.mat')==2
        disp('loading EXPT.mat'); load EXPT; 
    else
        disp('You need to create a variable that has information about the experiment (EXPT).');
        disp('This variable is called EXPT, and is stored in the main study directory in EXPT.mat')
        disp('No EXPT file in current directory.  Click Create to make one.');
        EXPT = [];
        fprintf(1,'\n')
    end
        
    varargout{1} = EXPT; 
        
        
    case lower('AsciiWelcome')
    %=======================================================================
    disp( 'Welcome to BrainVowager_gui.')
    fprintf('\n')

    case lower('Ver')
    %=======================================================================
    varargout = {'Tor Brainvowager V1'};
  


    case lower('CreateMenuWin')
    %=======================================================================
    close(findobj(get(0,'Children'),'Tag','BrainVowager_gui Menu'))

    
    % default sizes, etc.
    S = get(0,'ScreenSize');
    winwid = 300;               % window width
    winh = 400;                 % window height
    topbutton = winh-100;        % y location of top button
    butspace = 35;               % y spacing of buttons
    
    fullbutxy = [160 30];       % full-length button width and height
    rightbutxy = [100 30];      % right-hand button width and height
    halfbutxy = [60 30];        % (left-hand) half-width button w and h
    
    
    %-Open BrainVowager_gui menu window
    %-----------------------------------------------------------------------
    
    F = figure('Color',[1 1 1]*.8,...
        'Name',BrainVowager_gui('Ver'),...
        'NumberTitle','off',...
        'Position',[S(3)/2-450,S(4)/2-140,winwid,winh],...
        'Resize','off',...
        'Tag','BrainVowager_gui Menu',...
        'Pointer','Watch',...
        'MenuBar','none',...
        'Visible','off');

    %-Frames and text
    %-----------------------------------------------------------------------
    axes('Position',[0 0 80/winwid winh/winh],'Visible','Off')
    text(0.5,0.475,'BrainTools',...
        'FontName','Times','FontSize',36,...
        'Rotation',90,...
        'VerticalAlignment','middle','HorizontalAlignment','center',...
        'Color',[1 1 1]*.6);

    text(0.2,0.96,'SCN Lab',...
        'FontName','Times','FontSize',16,'FontAngle','Italic',...
        'FontWeight','Bold',...
        'Color',[1 1 1]*.6);

    uicontrol(F,'Style','Frame','Position',[095 005 winwid-100 winh - 30],...
        'BackgroundColor',BrainVowager_gui('Color'));  % colored frame
    uicontrol(F,'Style','Frame','Position',[105 015 winwid-120 winh - 50]);  % inner gray frame

    %-Buttons to launch BrainVowager_gui functions
    %-----------------------------------------------------------------------
    uicontrol(F,'Style','Text',...
        'String','Subject level','FontSize',14,...
        'HorizontalAlignment','Center',...
        'Position',[115 topbutton+30 fullbutxy],...
        'ForegroundColor','y','FontWeight','b');


    % FIR setup
    str = 'CreateExpt(''FIR'');';          % callback function
    str = BrainVowager_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Setup',...
        'Position',[110 topbutton-(1-1)*butspace 60 030],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    % FIR analysis
    str = 'EXPT = whole_brain_fir(EXPT,''full'',1);';       % callback function
    str = BrainVowager_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','FIR model',...
        'Position',[180 topbutton-(1-1)*butspace rightbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');

    % -------------------------------------------
        
    % EWMA setup
    str = 'CreateExpt(''hewma'');';                          % callback function
    str = BrainVowager_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Setup',...
        'Position',[110 topbutton-(2-1)*butspace 60 030],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    % EWMA analysis
    str = 'EXPT = wb_multisubject_ewma(EXPT);';             % callback function
    str = BrainVowager_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','EWMA',...
        'Position',[180 topbutton-(2-1)*butspace rightbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    % -------------------------------------------

    % Parametric setup
    str = 'CreateExpt(''parammod'');';                       % callback function
    str = BrainVowager_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Setup',...
        'Position',[110 topbutton-(3-1)*butspace 60 030],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    % Parametric analysis
    str = 'EXPT = wb_multisubject_parammod(EXPT);';         % callback function
    str = BrainVowager_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Parametric',...
        'Position',[180 topbutton-(3-1)*butspace rightbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    % -------------------------------------------
    % Next section
    % -------------------------------------------
    
    uicontrol(F,'Style','Text',...
    'String','Group level','FontSize',14,...
    'HorizontalAlignment','Center',...
    'Position',[115 topbutton-(4-1)*butspace fullbutxy],...
    'ForegroundColor','y','FontWeight','b');

    % -------------------------------------------

    % Robust setup
%     str = 'CreateExpt(''robfit'');';                       % callback function
%     str = BrainVowager_gui('ExpandString',str);                 % display then execute
%     uicontrol(F,'String','Setup',...
%         'Position',[110 topbutton-(5-1)*butspace 60 030],...
%         'CallBack',str,...
%         'Interruptible','on',...
%         'ForegroundColor','k','FontWeight','b');
    
    % Robust GUI
    str = 'scn_fir_results_gui;';                                % callback function
    str = BrainVowager_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','FIR RFX/Seed',...
        'Position',[110 topbutton-(5-1)*butspace fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');

    % -------------------------------------------
    % Next section
    % -------------------------------------------
    
    uicontrol(F,'Style','Text',...
    'String','ROI Tools','FontSize',14,...
    'HorizontalAlignment','Center',...
    'Position',[115 topbutton-(6-1)*butspace fullbutxy],...
    'ForegroundColor','y','FontWeight','b');

    % -------------------------------------------
    
    % ROI GUI
    str = 'scn_roi_gui;';                                % callback function
    str = BrainVowager_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','ROI Menu',...
        'Position',[180 topbutton-(7-1)*butspace rightbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    % -------------------------------------------



    set(F,'Pointer','Arrow','Visible','on')



    case lower('Color')
    %=======================================================================
    % BrainVowager_gui('Color')
    %-----------------------------------------------------------------------
    % %-Developmental livery
    % varargout = {[0.7,1.0,0.7], 'Lime Green'};
    %-Distribution livery
    varargout = {[.7 0.6 .5], 'Tan'};

    
    case lower('ExpandString')
    %=======================================================================
    % BrainVowager_gui('ExpandString')
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
        
        
    otherwise
    %=======================================================================
    error('Unknown action string')

    %======================================================================
    
end


return


