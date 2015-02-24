function varargout=scn_fir_results_gui(Action,varargin)
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
    %=======================================================================
    global EXPT
    %clc
    %BrainVowager_defaults;
    scn_fir_results_gui('AsciiWelcome')
    scn_fir_results_gui('CreateMenuWin')

    % load EXPT
    if ~(exist('EXPT.mat')==2)
        disp('You need to go to the directory containing individual subject results dirs ');
        disp('this directory should have information about the experiment stored in the file EXPT.mat')
        disp('No EXPT file in current directory.');
        EXPT = [];
        fprintf(1,'\n')
    end
    if ~isempty(EXPT)
        disp('Using EXPT already in memory.'); 
    elseif exist('EXPT.mat')==2
        disp('loading EXPT.mat'); load EXPT; 
    end
        
    varargout{1} = EXPT; 
        
        
        
    case lower('AsciiWelcome')
    %=======================================================================
    disp( 'Welcome to scn_fir_results: Contrasts and results for whole_brain_fir.')
    fprintf('\n')

    case lower('Ver')
    %=======================================================================
    varargout = {'SCNlab FIR results'};
  


    case lower('CreateMenuWin')
    %=======================================================================
    close(findobj(get(0,'Children'),'Tag','scn_fir_results_gui Menu'))

    %-Open scn_fir_results menu window
    %-----------------------------------------------------------------------
    [F, winwid, winh] = scn_fir_results_gui('initFigure');

    
    % default button sizes and positions, etc.
    
    topbutton = winh-100;        % y location of top button
    butspace = 35;               % y spacing of buttons
    
    fullbutxy = [160 30];       % full-length button width and height
    halfbutxy = [80 30];        % (left-hand) half-width button w and h
    rightbutx = 110+halfbutxy(1)+5;  % right-hand button start x
    
    
    %-Frames and text
    %-----------------------------------------------------------------------
    axes('Position',[0 0 80/winwid winh/winh],'Visible','Off')
    text(0.5,0.475,'FIR model results',...
        'FontName','Times','FontSize',36,...
        'Rotation',90,...
        'VerticalAlignment','middle','HorizontalAlignment','center',...
        'Color',[1 1 1]*.6);

    text(0.2,0.96,'Whole brain tools',...
        'FontName','Times','FontSize',16,'FontAngle','Italic',...
        'FontWeight','Bold',...
        'Color',[1 1 1]*.6);

    uicontrol(F,'Style','Frame','Position',[095 005 winwid-100 winh - 30],...
        'BackgroundColor',scn_fir_results_gui('Color'));  % colored frame
    uicontrol(F,'Style','Frame','Position',[105 015 winwid-120 winh - 50]);  % inner gray frame
    
    

    %-Buttons to launch scn_fir_results functions
    %-----------------------------------------------------------------------
    uicontrol(F,'Style','Text',...
        'String','Set up contrasts','FontSize',14,...
        'HorizontalAlignment','Center',...
        'Position',[115 topbutton+30 160 030],...
        'ForegroundColor','y','FontWeight','b');

    % Collect image names
    str = 'EXPT = get_htw_image_names(EXPT);';                 % callback function
    str = scn_fir_results_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Get H/T/W images',...
        'Position',[110 topbutton-(1-1)*35 160 030],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    

    % Specify contrasts
    str = 'EXPT = make_htw_contrast_images(EXPT);';            % callback function
    str = scn_fir_results_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Specify contrasts',...
        'Position',[110 topbutton-(2-1)*35 160 030],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    
    % -------------------------------------------
    
    uicontrol(F,'Style','Text',...
    'String','Group level','FontSize',14,...
    'HorizontalAlignment','Center',...
    'Position',[115 topbutton-(3-1)*35 160 030],...
    'ForegroundColor','y','FontWeight','b');

    % -------------------------------------------

    % Robfit
    %str = 'EXPT = robfit(EXPT, [], 0, EXPT.mask);';            % callback function
    str = 'robust_toolbox_gui;';
    str = scn_fir_results_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Robfit',...
        'Position',[110 topbutton-(4-1)*35 160 030],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');

    %Batch time-course
    %str = ['EXPT = CreateExpt(''CreateField'',''dxnames'',''FIR'');', ...   
    %    'EXPT = robust_batch_dx_plot(EXPT,''both'');' ];       % callback function

    str = 'EXPT = robust_batch_dx_plot(EXPT,''both'');';
    str = scn_fir_results_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String','Batch time-course',...
        'Position',[110 topbutton-(5-1)*35 160 030],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    % -------------------------------------------

    
    % -------------------------------------------

%     % Robust setup
%     str = 'CreateExpt(''robfit'');';                           % callback function
%     str = scn_fir_results('ExpandString',str);                 % display then execute
%     uicontrol(F,'String','Setup',...
%         'Position',[110 topbutton-(5-1)*35 60 030],...
%         'CallBack',str,...
%         'Interruptible','on',...
%         'ForegroundColor','k','FontWeight','b');
%     
%     % Parametric analysis
%     str = 'EXPT = robfit(EXPT);';                              % callback function
%     str = scn_fir_results_gui('ExpandString',str);                 % display then execute
%     uicontrol(F,'String','Robust RFX',...
%         'Position',[180 topbutton-(5-1)*35 100 030],...
%         'CallBack',str,...
%         'Interruptible','on',...
%         'ForegroundColor','k','FontWeight','b');
    
    % -------------------------------------------



    set(F,'Pointer','Arrow','Visible','on')



    case lower('Color')
    %=======================================================================
    % scn_fir_results_gui('Color')
    %-----------------------------------------------------------------------
    % 
    
    varargout = {[1 0.8 0.6], 'Orange'};

    
    case lower('ExpandString')
    %=======================================================================
    % scn_fir_results_gui('ExpandString')
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
    % [F, winwid, winh] = scn_fir_results_gui('initFigure')
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
    
    if ~isempty(h), 
        pos = get(h,'Position');, 
        winwid = pos(3); winh = pos(4);
        pos(1) = pos(1) + winwid;   % put next to main figure
    end
    
    %-Open scn_fir_results_gui menu window
    %----------------------------------------------------------------------
    
    F = figure('Color',[1 1 1]*.8,...
        'Name',scn_fir_results_gui('Ver'),...
        'NumberTitle','off',...
        'Position',pos,...
        'Resize','off',...
        'Tag','scn_fir_results_gui Menu',...
        'Pointer','Watch',...
        'MenuBar','none',...
        'Visible','off');
    
    varargout{1} = F; varargout{2} = winwid; varargout{3} = winh;
    
    
     
        
    otherwise
    %=======================================================================
    error('Unknown action string')

    %======================================================================
    
end


return


