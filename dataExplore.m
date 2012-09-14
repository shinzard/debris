function varargout = dataExplore(varargin)
% DATAEXPLORE M-file for dataExplore.fig
%      DATAEXPLORE, by itself, creates a new DATAEXPLORE or raises the existing
%      singleton*.
%
%      H = DATAEXPLORE returns the handle to a new DATAEXPLORE or the handle to
%      the existing singleton*.
%
%      DATAEXPLORE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATAEXPLORE.M with the given input arguments.
%
%      DATAEXPLORE('Property','Value',...) creates a new DATAEXPLORE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dataExplore_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dataExplore_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dataExplore

% Last Modified by GUIDE v2.5 30-Aug-2011 15:08:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dataExplore_OpeningFcn, ...
                   'gui_OutputFcn',  @dataExplore_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before dataExplore is made visible.
function dataExplore_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dataExplore (see VARARGIN)

% Choose default command line output for dataExplore
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using dataExplore.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes dataExplore wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dataExplore_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in updatePlot.
function updatePlot_Callback(hObject, eventdata, handles)
% hObject    handle to updatePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.dataAx);
cla;

data_index = get(handles.dataType, 'Value');
category_index = get(handles.category, 'Value');
select_indx = get(handles.select, 'Value');
switch data_index
  case 1 %GPS
    AL = shaperead('/home/james/Documents/Research/Debris Task Assignment/GIS/ALcounties/tl_2010_01_cousub10.shp', 'UseGeoCoords', true);
    hold off;
    for i = 1:length(AL)
        if sum(AL(i).BoundingBox(:,2) > 32) > 1
            plot(AL(i).Lon, AL(i).Lat, 'k--', 'LineWidth', 0.5); hold on;    
        end
    end

    axis tight;
    axis equal;
    axis manual;
    %axis off;
    
    lat = evalin('base', 'lat');
    lon = evalin('base', 'lon');
    
    switch category_index
      case 1 % Dummy
      case 2 % TruckID
        truckId = evalin('base', 'truckId');
        trucks = evalin('base', 'trucks');  
        idx = find(lon ~= 0);
        idx2 = (truckId(idx) == trucks(select_indx));
        if get(handles.includeOtherData, 'Value')
            plot(lon(idx), lat(idx), 'b.', 'MarkerSize', 4);
            hold on;
        end
        
        plot(lon(idx(idx2)), lat(idx(idx2)), 'r.')
        title(sprintf('Ticket data for truck %d', trucks(select_indx)), 'FontSize', 12);
        
      case 3 % QA
        QA = evalin('base', 'QA');
        QAs = unique(QA);
        idx = find(lon ~= 0);
        idx2 = (QA(idx) == QAs(select_indx));
        if get(handles.includeOtherData, 'Value')
            plot(lon(idx), lat(idx), 'b.', 'MarkerSize', 4);
            hold on;
        end
        
        plot(lon(idx(idx2)), lat(idx(idx2)), 'r.')
        title(sprintf('Ticket data for QA %d', QAs(select_indx)), 'FontSize', 12);
        
      case 4 % QC
        QC = evalin('base', 'QC');
        QCs = unique(QC);
        idx = find(lon ~= 0);
        idx2 = (QC(idx) == QCs(select_indx));
        if get(handles.includeOtherData, 'Value')
            plot(lon(idx), lat(idx), 'b.', 'MarkerSize', 4);
            hold on;
        end
        
        plot(lon(idx(idx2)), lat(idx(idx2)), 'r.')
        title(sprintf('Ticket data for QC %d', QCs(select_indx)), 'FontSize', 12);
      case 5 % Sub
        subcont = evalin('base', 'subcont');
        subconts = unique(subcont);
        idx = find(lon ~= 0);
        idx2 = (subcont(idx) == subconts(select_indx));
        if get(handles.includeOtherData, 'Value')
            plot(lon(idx), lat(idx), 'b.', 'MarkerSize', 4);
            hold on;
        end
        
        plot(lon(idx(idx2)), lat(idx(idx2)), 'r.')
        title(sprintf('Ticket data for Subcontractor %d', subconts(select_indx)), 'FontSize', 12);
                
      case 6 %Region
        tdsr = evalin('base', 'tdsr');
        tdsrs = unique(tdsr);
        
        idx = find(lon ~= 0);
        idx2 = (tdsr(idx) == tdsrs(select_indx));
        if get(handles.includeOtherData, 'Value')
            plot(lon(idx), lat(idx), 'b.', 'MarkerSize', 4);
            hold on;
        end
        
        plot(lon(idx(idx2)), lat(idx(idx2)), 'r.')
        title(sprintf('Ticket data for TDSR %d', ...
                      tdsrs(select_indx)), 'FontSize', 12);
      case 7  %Project
        project = evalin('base', 'project');
        projects = unique(project);
        idx = find(lon ~= 0);
        idx2 = (project(idx) == projects(select_indx));
        if get(handles.includeOtherData, 'Value')
            plot(lon(idx), lat(idx), 'b.', 'MarkerSize', 4);
            hold on;
        end
        
        plot(lon(idx(idx2)), lat(idx(idx2)), 'r.')
        title(sprintf('Ticket data for project (municipality) %2.1f', ...
                      projects(select_indx)), 'FontSize', 12);
      case 8  %Debris Type
        contents = evalin('base', 'contents');
        type = unique(contents);
        idx = find(lon ~= 0);
        idx2 = (contents(idx) == type(select_indx));
        if get(handles.includeOtherData, 'Value')
            plot(lon(idx), lat(idx), 'b.', 'MarkerSize', 4);
            hold on;
        end
        
        plot(lon(idx(idx2)), lat(idx(idx2)), 'r.')
        title(sprintf('Ticket data for debris type %2.1f', ...
                      type(select_indx)), 'FontSize', 12);
        
    end
        
    case 2 % Other data...

end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in dataType.
function dataType_Callback(hObject, eventdata, handles)
% hObject    handle to dataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns dataType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dataType


% --- Executes during object creation, after setting all properties.
function dataType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'GPS', 'LoadData'});


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over updatePlot.
function updatePlot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to updatePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in category.
function category_Callback(hObject, eventdata, handles)
% hObject    handle to category (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns category contents as cell array
%        contents{get(hObject,'Value')} returns selected item from category

category_index = get(handles.category, 'Value');

switch category_index
  case 1 % Dummy
  case 2 % Trucks
    trucks = evalin('base', 'trucks');
    set(handles.select, 'String', num2cell(trucks));
  case 3 % QA
    QA = evalin('base', 'QA');
    set(handles.select, 'String', num2cell(unique(QA)));
  case 4 % QC
    QC = evalin('base', 'QC');
    set(handles.select, 'String', num2cell(unique(QC)));
  case 5 % Sub
    subcont = evalin('base', 'subcont');
    set(handles.select, 'String', num2cell(unique(subcont)));
  case 6 % TDSR
    tdsr = evalin('base', 'tdsr');
    set(handles.select, 'String', num2cell(unique(tdsr)));
  case 7 % Project
    project = evalin('base', 'project');
    set(handles.select, 'String', num2cell(unique(project)));
  case 8 % Debris type
    contents = evalin('base', 'contents');
    num2cell(unique(contents))
    set(handles.select, 'String', num2cell(unique(contents)));
end


% --- Executes during object creation, after setting all properties.
function category_CreateFcn(hObject, eventdata, handles)
% hObject    handle to category (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in select.
function select_Callback(hObject, eventdata, handles)
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select


% --- Executes during object creation, after setting all properties.
function select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in includeOtherData.
function includeOtherData_Callback(hObject, eventdata, handles)
% hObject    handle to includeOtherData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of includeOtherData
