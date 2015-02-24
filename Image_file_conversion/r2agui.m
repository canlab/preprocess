function varargout = r2agui(varargin)
% R2AGUI Application M-file for r2agui.fig
%
%   This program converts the Philips PAR and REC image file format
%   to the Analyze or Nifti 1.0 format used by SPM5 and other image processing packages
%   Works for the very recent PAR & REC file version V3 and V4 
%   Type 'r2agui' on the matlab command line to invoke the graphical user
%   interface. Make sure that r2agui.m and r2agui.fig as well as the SPM5 functions are in the
%   matlab search path.
%
%   Created by:
%   Erno Hermans and Bas Neggers, Helmholtz Institute,
%   Utrecht University/Dept of Brain Research, University Medical Center
%
%   Version 2.2, released 31-10-2006
global r2a_defaults;


if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename,'reuse');

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    guidata(fig, handles);
    
    %load r2agui default settings, when found
    if exist('r2agui_ini.m')
        r2agui_ini;
    else
        disp('no r2agui_ini with user specific settings found,using r2agui_defaults');
        r2agui_defaults;
    end
    
    %initialize some handles according to defaults
    if exist(r2a_defaults.altfolder)
        set(handles.usealtfolder, 'Value',r2a_defaults.usealtfolder);
        set(handles.altfolder,'String',r2a_defaults.altfolder);
    end
    set(handles.outputformat,'value',r2a_defaults.outputformat);
    
    
    if nargout > 0
        varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and
%| sets objects' callback properties to call them through the FEVAL
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = batchtoggle_Callback(h, eventdata, handles, varargin)

switch(get(handles.batchtoggle,'value')),
    case 1
        disp('Process all volumes in directory');
        set(handles.batchtoggle,'string','Batch on');
    case 0
        disp('Process selected volume only');
        set(handles.batchtoggle,'string','Batch off');
end

% --------------------------------------------------------------------

function varargout = openpar_Callback(h, eventdata, handles, varargin)
global pathpar
global file
global r2a_defaults
pathpar=r2a_defaults.PARREC_path;


parfile = r2aselect(1,pathpar,'select PAR file','.PAR');
[pathpar, filename,ext]=fileparts(parfile);
pathpar=[pathpar,filesep];
file=[filename,ext];
Pars=read_par(parfile);

set(handles.voxelsizex, 'String',Pars.vox(1));
set(handles.voxelsizey, 'String',Pars.vox(2));
set(handles.voxelsizez, 'String',Pars.vox(3));
set(handles.dynamics, 'String',Pars.dyn);
set(handles.scanduration, 'String',Pars.RT);
set(handles.RT_Version, 'String',Pars.ResToolsVersion);
set(handles.filenaam, 'String',parfile);
prefix = get(handles.prefix, 'String');

subaan = get(handles.subdir, 'Value');
if subaan ==1
    if get(handles.fullprefix,'value')==1
        subdir=prefix;
    else
        subdir = [prefix,file(1:(length(file)-4))];
    end
else
    subdir = '';
end
setfoldernames(handles);



set(handles.convert, 'Enable', 'on');

% --------------------------------------------------------------------
function varargout = voxelsizex_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = voxelsizey_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = voxelsizez_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = dynamics_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = filenaam_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = prefix_Callback(h, eventdata, handles, varargin)
setfoldernames(handles);

% --------------------------------------------------------------------
function varargout = subdir_Callback(h, eventdata, handles, varargin)
setfoldernames(handles);

% --------------------------------------------------------------------
function varargout = usealtfolder_Callback(h, eventdata, handles, varargin)
setfoldernames(handles);

% --------------------------------------------------------------------
function varargout = altfolder_Callback(h, eventdata, handles, varargin)
setfoldernames(handles);

function varargout=setfoldernames(handles)
global pathpar
global file
if ispc,
    WinLinSlash='\';
else
    WinLinSlash='/';
end

niioranalyze={'nii','img'};
filext=niioranalyze{get(handles.outputformat,'value')};

prefix = get(handles.prefix, 'String');
subaan = get(handles.subdir, 'Value');
altfolderaan=get(handles.usealtfolder, 'Value');


if altfolderaan ==1;
    outfoldername = get(handles.altfolder, 'String');
    i=find(pathpar==filesep);
    if isempty(outfoldername)
        disp('No alternate folder specified, using current + name lowest level subdirectory containing parfile');
        outfoldername=pwd;
    end
    if ~strcmp(outfoldername(end),WinLinSlash) %check if foldername has a trailing slash, add when necessary
        outfoldername=[outfoldername,WinLinSlash];
    end
    set(handles.altfolder, 'String',outfoldername); % dont show extra PAR file folder here, but add to path below
    outfoldername=[outfoldername,pathpar(i(end-1)+1:i(end))];
else
    outfoldername = pathpar;
end

if subaan ==1
    if get(handles.fullprefix,'value')==1
        subdir=prefix;
    else
        subdir = [prefix,file(1:(length(file)-4))];
    end
else
    subdir = '';
end
if get(handles.fullprefix,'value')==1
    set(handles.outputfile, 'String',fullfile(outfoldername,subdir, [prefix,'-001.',filext]));
else
    set(handles.outputfile, 'String',fullfile(outfoldername,subdir, [prefix,file(1:(length(file)-4)),'-001.',filext]));
end


% --------------------------------------------------------------------
function varargout = edit8_Callback(h, eventdata, handles, varargin)







% --------------------------------------------------------------------r2a

function varargout = convert_Callback(h, eventdata, handles, varargin)
global file;
global pathpar;
global r2a_defaults;

if ispc,
    WinLinSlash='\';
else
    WinLinSlash='/';
end
if get(handles.batchtoggle,'value')==0
    no_files = 1;
    filelist{1}=file;
else
    if ispc,    %check whether this is windows
        dl=dir([pathpar,'*.par']);
        no_files = length(dl);
        for i=1:length(dl),
            filelist{i}=dl(i).name;
        end
    else
        dlL=dir([pathpar,'*.par']);
        dlU=dir([pathpar,'*.PAR']);
        no_files = length(dlL)+length(dlU);
        for i=1:length(dlL),
            filelist{i}=dlL(i).name;
        end
        if length(dlL)==0,
            i=0;
        end
        for i2=1:length(dlU),
            filelist{i+i2}=dlU(i2).name;
        end
    end
end

options.subaan=get(handles.subdir,'Value');
options.usealtfolder=get(handles.usealtfolder, 'Value');
options.altfolder=get(handles.altfolder, 'String');
options.prefix=get(handles.prefix, 'String');
options.pathpar=pathpar;
options.angulation=r2a_defaults.angulation;
options.rescale=r2a_defaults.rescale;
options.usefullprefix=get(handles.fullprefix,'value');
options.outputformat=get(handles.outputformat,'value');
outfile=convert_r2a(filelist,options);

set(handles.outputfile, 'String',outfile{1});
set(handles.convert, 'Enable', 'on');


% --------------------------------------------------------------------
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)






function RT_Version_Callback(hObject, eventdata, handles)
% hObject    handle to RT_Version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RT_Version as text
%        str2double(get(hObject,'String')) returns contents of RT_Version as a double


% --- Executes during object creation, after setting all properties.
function RT_Version_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RT_Version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function altfolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to altfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in outputformat.
function outputformat_Callback(hObject, eventdata, handles)
% hObject    handle to outputformat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns outputformat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from outputformat


% --- Executes during object creation, after setting all properties.
function outputformat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputformat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fullprefix.
function fullprefix_Callback(hObject, eventdata, handles)
% hObject    handle to fullprefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fullprefix
setfoldernames(handles);


function outfiles=r2aselect(nr,spath,stitle,swild);
%spm-version aware file selector, because in X systems default matlab file
%selector sucks

if exist('spm') % check of any spm in searchpatch
    switch(spm('ver')) % ask spm which version it is
        case 'SPM2'
           outfiles=spm_get(nr,swild,stitle,spath);
        case 'SPM5'
           outfiles=spm_select(nr,'any',stitle,'',spath,swild,1);
        otherwise
           display([spm('ver') ' version of SPM software was unknown at the time this r2agui was created']);
    end
else % no spm anywhere near (r2agui started in some scary place)
    disp('no SPM version found, using default matlab file selector (yikes!)');
    [file,pathfile]=uigetfile(['*' swild],stitle,'MultiSelect', 'on'); 
    outfiles=[pathfile,file];
end



function scanduration_Callback(hObject, eventdata, handles)
% hObject    handle to scanduration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanduration as text
%        str2double(get(hObject,'String')) returns contents of scanduration as a double


% --- Executes during object creation, after setting all properties.
function scanduration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scanduration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


