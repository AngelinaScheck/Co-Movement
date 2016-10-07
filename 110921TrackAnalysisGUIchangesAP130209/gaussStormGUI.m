
function varargout = gaussStormGUI(varargin)
% GAUSSSTORMGUI MATLAB code for gaussStormGUI.fig
%      GAUSSSTORMGUI, by itself, creates a new GAUSSSTORMGUI or raises the existing
%      singleton*.
%
%      H = GAUSSSTORMGUI returns the handle to a new GAUSSSTORMGUI or the handle to
%      the existing singleton*.
%
%      GAUSSSTORMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAUSSSTORMGUI.M with the given input arguments.
%
%      GAUSSSTORMGUI('Property','Value',...) creates a new GAUSSSTORMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gaussStormGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gaussStormGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gaussStormGUI

% Last Modified by GUIDE v2.5 07-Oct-2016 13:21:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gaussStormGUI_OpeningFcn, ...
    'gui_OutputFcn',  @gaussStormGUI_OutputFcn, ...
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


% --- Executes just before gaussStormGUI is made visible.
function gaussStormGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gaussStormGUI (see VARARGIN)

% Choose default command line output for gaussStormGUI
handles.output = hObject;
% use guiMain to initialise appData
guiMain('init',handles, varargin);
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = gaussStormGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%% create functions

function trackingWindow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotTracksMinSteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotTracksColour_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pixel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sigmaNoise_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bacLength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bacWidth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxStep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rangeD_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DhistMinSteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function initGuess_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CDFfittingOptions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function localizationThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function greenThreshold_CreateFcn(hObject, eventdata, handles) %AS for co-movement
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function redThreshold_CreateFcn(hObject, eventdata, handles) %AS for co-movement
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DhistThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function curveColor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FilterLnoise_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FilterLobject_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackMemory_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numberOfSubplots_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TracksSelection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PlotTracksSelection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function findFilterParamsOption_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minImmobilizedFrames_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filterDwindow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function AltType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function startframe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function leftRim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%% callbacks

function tracking_Callback(hObject, eventdata, handles)
guiMain('tracking_Callback',handles);

function trackingWindow_Callback(hObject, eventdata, handles)
guiMain('trackingWindow_Callback',handles);

function plotTracks_Callback(hObject, eventdata, handles)
guiMain('plotTracks_Callback',handles);

function plotTracksMinSteps_Callback(hObject, eventdata, handles)
guiMain('plotTracksMinSteps_Callback',handles);

function plotTracksColour_Callback(hObject, eventdata, handles)
guiMain('plotTracksColour_Callback',handles);

function checkboxSaveTracks_Callback(hObject, eventdata, handles)
guiMain('checkboxSaveTracks_Callback',handles);

function MSDcurve_Callback(hObject, eventdata, handles)
guiMain('MSDcurve_Callback',handles);

function pixel_Callback(hObject, eventdata, handles)
guiMain('pixel_Callback',handles);

function dT_Callback(hObject, eventdata, handles)
guiMain('dT_Callback',handles);

function sigmaNoise_Callback(hObject, eventdata, handles)
guiMain('sigmaNoise_Callback',handles);

function bacLength_Callback(hObject, eventdata, handles)
guiMain('bacLength_Callback',handles);

function bacWidth_Callback(hObject, eventdata, handles)
guiMain('bacWidth_Callback',handles);

function maxStep_Callback(hObject, eventdata, handles)
guiMain('maxStep_Callback',handles);

function clearData_Callback(hObject, eventdata, handles)
guiMain('clearData_Callback',handles);

function histD_Callback(hObject, eventdata, handles)
guiMain('histD_Callback',handles);

function rangeD_Callback(hObject, eventdata, handles)
guiMain('rangeD_Callback',handles);

function DhistMinSteps_Callback(hObject, eventdata, handles)
guiMain('DhistMinSteps_Callback',handles);

function CDFcurve_Callback(hObject, eventdata, handles)
guiMain('CDFcurve_Callback',handles);

function initGuess_Callback(hObject, eventdata, handles)
guiMain('initGuess_Callback',handles);

function CDFfittingOptions_Callback(hObject, eventdata, handles)
guiMain('CDFfittingOptions_Callback',handles);

function plotLocalizations_Callback(hObject, eventdata, handles)
guiMain('plotLocalizations_Callback',handles);

function selectROIs_Callback(hObject, eventdata, handles)
guiMain('selectROIs_Callback',handles);

function localization_Callback(hObject, eventdata, handles)
guiMain('localization_Callback',handles);

function localizationThresh_Callback(hObject, eventdata, handles)
guiMain('localizationThresh_Callback',handles);

function findBindingDiffusionTracks_Callback(hObject, eventdata, handles)
guiMain('findBindingDiffusionTracks_Callback',handles);

function DThresh_Callback(hObject, eventdata, handles)
guiMain('DThresh_Callback',handles);

function TwoSpeciesMSDThreshold_Callback(hObject, eventdata, handles)
guiMain('TwoSpeciesMSDThreshold_Callback',handles);

function DhistThresh_Callback(hObject, eventdata, handles)
guiMain('DhistThresh_Callback',handles);

function holdFigureCheckbox_Callback(hObject, eventdata, handles)
guiMain('holdFigureCheckbox_Callback',handles);

function curveColor_Callback(hObject, eventdata, handles)
guiMain('curveColor_Callback',handles);

function combineTracks_Callback(hObject, eventdata, handles)
guiMain('combineTracks_Callback',handles);

function plotBrightfield_Callback(hObject, eventdata, handles)
guiMain('plotBrightfield_Callback',handles);

function checkboxAviMovie_Callback(hObject, eventdata, handles)
guiMain('checkboxAviMovie_Callback',handles);

function checkboxTifStack_Callback(hObject, eventdata, handles)
guiMain('checkboxTifStack_Callback',handles);

function checkboxFindFilterParams_Callback(hObject, eventdata, handles)
guiMain('checkboxFindFilterParams_Callback',handles);

function checkboxLoadTifFile_Callback(hObject, eventdata, handles)
guiMain('checkboxLoadTifFile_Callback',handles);

function FilterLnoise_Callback(hObject, eventdata, handles)
guiMain('FilterLnoise_Callback',handles);

function FilterLobject_Callback(hObject, eventdata, handles)
guiMain('FilterLobject_Callback',handles);

function FRETanalysis_Callback(hObject, eventdata, handles)
guiMain('FRETanalysis_Callback',handles);

function checkboxfirstframeG_Callback(hObject, eventdata, handles)
guiMain('checkboxfirstframeG_Callback',handles);

function trackMemory_Callback(hObject, eventdata, handles)
guiMain('trackMemory_Callback',handles)

function numberOfSubplots_Callback(hObject, eventdata, handles)
guiMain('numberOfSubplots_Callback',handles)

function TracksSelection_Callback(hObject, eventdata, handles)
guiMain('TracksSelection_Callback',handles)

function PlotTracksSelection_Callback(hObject, eventdata, handles)
guiMain('PlotTracksSelection_Callback',handles)

function findFilterParamsOption_Callback(hObject, eventdata, handles)
guiMain('findFilterParamsOption_Callback',handles)

function UBTimeTrace_Callback(hObject, eventdata, handles)
guiMain('UBTimeTrace_Callback',handles)

function filterDwindow_Callback(hObject, eventdata, handles)
guiMain('filterDwindow_Callback',handles)

function minImmobilizedFrames_Callback(hObject, eventdata, handles)
guiMain('minImmobilizedFrames_Callback',handles)

function checkboxDallSelect_Callback(hObject, eventdata, handles)
guiMain('checkboxDallSelect_Callback',handles)

function AltType_Callback(hObject, eventdata, handles)
guiMain('AltType_Callback',handles)

function startframe_Callback(hObject, eventdata, handles)
guiMain('startframe_Callback',handles)

function plotFRETdata_Callback(hObject, eventdata, handles)
guiMain('plotFRETdata_Callback',handles)

function overlaymovie_Callback(hObject, eventdata, handles)
guiMain('overlaymovie_Callback',handles)

function FRETROIs_Callback(hObject, eventdata, handles)
guiMain('FRETROIs_Callback',handles)

function FRETtracking_Callback(hObject, eventdata, handles)
guiMain('FRETtracking_Callback',handles)

function leftRim_Callback(hObject, eventdata, handles)
guiMain('leftRim_Callback',handles)

%Added by AS for comovement analysis
function comovement_Callback(hObject, eventdata, handles)
guiMain('comovement_Callback', handles);

function redThreshold_Callback(hObject, eventdata, handles) %AS for co-movement
guiMain('redThreshold_Callback',handles);

function greenThreshold_Callback(hObject, eventdata, handles)%AS for co-movement
guiMain('greenThreshold_Callback',handles);

function comovTrackConnect_Callback(hObject, eventdata, handles)
guiMain('comovTrackConnect_Callback',handles);


%------------------------------------------------------------------------
%--------application specific code
function guiMain(param, handles,varargin)
% function guiMain(param, handles)
% ----------Main control function ---------------------------------------

if strcmp(param, 'init') % on initialise, create the appData variable & initialise fields
    
    appData.trackParams.mem = 0;
    appData.trackParams.dim = 2;
    appData.trackParams.good = 0;
    appData.trackParams.quiet = 0;
    appData.trackParams.maxDisp = 5;
    
    set(handles.trackingWindow,'String',num2str(appData.trackParams.maxDisp));
    set(handles.trackMemory,'String',num2str(appData.trackParams.mem));
    
    appData.plotTracksMinSteps = 4;
    set(handles.plotTracksMinSteps,'String',num2str(appData.plotTracksMinSteps));
    
    appData.numberOfSubplots = 6;
    set(handles.numberOfSubplots,'String',num2str(appData.numberOfSubplots));
    
    appData.plotTracksColour = 'time';
    set(handles.plotTracksColour,'String',appData.plotTracksColour);
    
    setappdata(handles.figure1, 'appData', appData);
    
    appData.checkboxSaveTracks = 1;
    set(handles.checkboxSaveTracks,'Value',appData.checkboxSaveTracks);
    
    appData.checkboxfirstframeG = 1;
    set(handles.checkboxfirstframeG,'Value',appData.checkboxfirstframeG);
        
    appData.checkboxAviMovie = 0;
    set(handles.checkboxAviMovie,'Value',appData.checkboxAviMovie);
    
    appData.checkboxTifStack = 0;
    set(handles.checkboxTifStack,'Value',appData.checkboxTifStack);
    
    appData.checkboxFindFilterParams = 0;
    set(handles.checkboxFindFilterParams,'Value',appData.checkboxFindFilterParams);
    
    appData.checkboxLoadTifFile = 0;
    set(handles.checkboxLoadTifFile,'Value',appData.checkboxLoadTifFile);
     
    appData.checkboxDallSelect = 0;
    set(handles.checkboxDallSelect,'Value',appData.checkboxDallSelect);
    
    appData.FilterLnoise = 1;
    set(handles.FilterLnoise,'String',num2str(appData.FilterLnoise));
    
    appData.FilterLobject = 7;
    set(handles.FilterLobject,'String',num2str(appData.FilterLobject));
    
    appData.pixel = 0.094;
    set(handles.pixel,'String',num2str(appData.pixel));
    
    appData.dT = 0.01028;
    set(handles.dT,'String',num2str(appData.dT));
    
    appData.sigmaNoise = 0.03/appData.pixel;
    set(handles.sigmaNoise,'String',num2str(appData.sigmaNoise*appData.pixel));
    
    appData.bacLength = 0.6; % length of rectangular section, not including round caps
    set(handles.bacLength,'String',num2str(appData.bacLength));
    
    appData.bacWidth = 0.6;
    set(handles.bacWidth,'String',num2str(appData.bacWidth));
    
    appData.bacArea = appData.bacLength * appData.bacWidth / (appData.pixel^2);
    
    appData.maxStep = 8;
    set(handles.maxStep,'String',num2str(appData.maxStep));
    
    appData.rangeDString = '-0.2:0.1:4'; % D range histD
    appData.rangeD = str2num(appData.rangeDString);
    set(handles.rangeD,'String',appData.rangeDString);
    
    appData.DhistMinSteps = 4;
    set(handles.DhistMinSteps,'String',num2str(appData.DhistMinSteps));
    
    appData.initGuess = 1;
    set(handles.initGuess,'String',num2str(appData.initGuess));
    
    appData.CDFfittingOption = 2;
    set(handles.CDFfittingOptions,'Value',appData.CDFfittingOption);
    
    appData.TracksSelection = 2;
    set(handles.TracksSelection,'Value',appData.TracksSelection);
    
    appData.findFilterParamsOption = 2;
    set(handles.findFilterParamsOption,'Value',appData.findFilterParamsOption);
    
    appData.PlotTracksSelection = 1;
    set(handles.TracksSelection,'Value',appData.TracksSelection);
    
    %appData.localizationWindow = 7; was set to defined value previously
    appData.localizationThresh = 20;
    set(handles.localizationThresh,'String',num2str(appData.localizationThresh));
    
    appData.redThreshold = 10;
    set(handles.redThreshold,'String',num2str(appData.redThreshold)); %AS for co-movement
    
    appData.greenThreshold = 10;
    set(handles.greenThreshold,'String',num2str(appData.greenThreshold)); %AS for co-movement
    
    appData.nFiles = 1;
    
    appData.DThreshString = '0.1 0.4'; %findBindingDiffusionTracks Threshold
    appData.DThresh = str2num(appData.DThreshString);
    set(handles.DThresh,'String',appData.DThreshString);
    
    appData.DhistThresh = 0.1; %TwoSpeciesMSDThreshold Threshold
    set(handles.DhistThresh,'String',num2str(appData.DhistThresh));
    
    appData.holdFigureCheckbox = 0;
    set(handles.holdFigureCheckbox,'Value',appData.holdFigureCheckbox);
    
    appData.curveColor = 'b';
    set(handles.curveColor,'String',appData.curveColor);
    
    appData.filterDwindow = 2; %filter D window
    set(handles.filterDwindow,'String',num2str(appData.filterDwindow));
    
    appData.minImmobilizedFrames = 8; %number of frames where mol is immobilized
    set(handles.minImmobilizedFrames,'String',num2str(appData.minImmobilizedFrames));
    
    appData.startframe = 1; %number of frames where mol is immobilized
    set(handles.startframe,'String',num2str(appData.startframe));
    
    appData.AltType = 1;
    set(handles.AltType,'Value',appData.AltType);
    
    appData.leftRim = 25;
    set(handles.leftRim,'String',num2str(appData.leftRim));
    
    
else  % subsequent runs, retrieve appData
    if isappdata(handles.figure1, 'appData')
        appData = getappdata(handles.figure1, 'appData');
    else
        error('appData, main data variable structure for avgGui not initialised');
    end
end


switch param
    
    case 'localization_Callback'
        
        [movieFilename, moviePathname] =...
            uigetfile({'*.fits';'*.tif'}, 'Movie data:','MultiSelect', 'on');
        
        if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
        
        info = whos('movieFilename');
        if strcmp(info.class,'char')
            appData.nFiles = 1;
        else
            appData.nFiles = numel(movieFilename);
        end
        
        for ii = 1:appData.nFiles
            
            if appData.nFiles == 1
                loadname  = [moviePathname movieFilename];
            else
                loadname = [moviePathname movieFilename{1,ii}];
            end
           
            if appData.checkboxFindFilterParams == 1
                FilterAverParam = findFilterParameters(loadname,appData.checkboxLoadTifFile,appData.findFilterParamsOption);
                appData.localizationThresh = FilterAverParam(1,1);
                appData.FilterLnoise = FilterAverParam(1,2);
                appData.FilterLobject = FilterAverParam(1,3);
            end
            
            gaussStormModifiedAP(loadname,appData.checkboxLoadTifFile,...
                appData.localizationThresh, appData.FilterLnoise, appData.FilterLobject);    
                        
        end
        
        if appData.nFiles == 1
            appData.data = importdata([moviePathname movieFilename(1:end-5),...
                '_thresh',num2str(appData.localizationThresh),...
                '_lnoise',num2str(appData.FilterLnoise),...
                '_lobject',num2str(appData.FilterLobject),'.out']);%load position data in workspace
        end
        appData.dataPathname = moviePathname;
        appData.dataFilename = movieFilename;
        end
            
        
    case 'localizationThresh_Callback'
        
        localizationThresh = str2double(get(handles.localizationThresh,'String'));
        
        if isnan(localizationThresh)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.localizationThresh = localizationThresh;
        
    case 'redThreshold_Callback' %AS for Co-Movement analysis
        
        RedThresh = str2num(get(handles.redThreshold,'String'));
        
        if isnan(RedThresh)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.redThreshold = RedThresh;
        
    case 'greenThreshold_Callback' %AS for Co-Movement analysis
        
        GreenThresh = str2num(get(handles.greenThreshold,'String'));
        
        if isnan(GreenThresh)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.greenThreshold = GreenThresh;
        
        
    case 'tracking_Callback'
        
        [dataFilename, dataPathname] =...
            uigetfile('*.out', 'Localization data:','MultiSelect', 'on');
        
        if ~(isnumeric(dataFilename)&&dataFilename==0) %check the user has not pressed cancel
        
        info = whos('dataFilename');
        if strcmp(info.class,'char')
            appData.nFiles = 1;
        else
            appData.nFiles = numel(dataFilename);
        end
        
        for ii = 1:appData.nFiles
            
            if appData.nFiles == 1
                loadname  = [dataPathname dataFilename];
                appData.dataPathname = dataPathname;
                appData.dataFilename = dataFilename;
            else
                loadname = [dataPathname dataFilename{1,ii}];
            end
            
            newData = importdata(loadname);
            if isstruct(newData)
                appData.data = newData.data;
            else
                appData.data = newData;
            end
            clear newData;
            
            %sort out positions depending on the 
            
            pos = zeros(length(appData.data(:,1)),3);
            
            % standard indexing
            pos(:,1) = appData.data(:,2);
            pos(:,2) = appData.data(:,3);
            pos(:,3) = appData.data(:,1);
            
            % old ROI indexing
            %             pos(:,1) = appData.data(:,1);
            %             pos(:,2) = appData.data(:,2);
            %             pos(:,3) = appData.data(:,4);
            
            tracks = trackWithDummy(pos, appData.trackParams);
            nMolecules = max(tracks(:,4));
            appData.tracks = tracks;
            appData.nMolecules = nMolecules;
            
            set(handles.results,'String',['nMolecules = ' num2str(nMolecules)]);
            
        
             if appData.checkboxSaveTracks == 1
%           choose where to save the file                     
%                 [tracksFilename,tracksPathname]=uiputfile([loadname '.tracks'],'Save tracking file as');
%                 loadnametracks = [tracksPathname tracksFilename];
%                 save([loadnametracks , '_win',num2str(appData.trackParams.maxDisp),'.tracks']);

%           save file in same directory with useful info like tracking window size
                save([loadname(1:end-4),'_win',num2str(appData.trackParams.maxDisp),'_mem',num2str(appData.trackParams.mem),'.tracks']);   
             end
            
            if appData.checkboxAviMovie == 1
            avimovie(tracks,appData.checkboxLoadTifFile,appData.plotTracksMinSteps);
            end;
            
            if appData.checkboxTifStack == 1
            tiffstack(tracks,appData.checkboxLoadTifFile,appData.plotTracksMinSteps);
            end

        end
        
        end
        
    case 'FRETanalysis_Callback'    %calls FRETSTORMAPModified
        
            %load movie data
        [movieFilename, moviePathname] =...
            uigetfile({'*.fits';'*.tif'}, 'Image data:','MultiSelect', 'on');
        
        if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
        
            %load Tform
        [TFormFilename, TFormPathname] =...
            uigetfile('*.mat', 'Transform .mat-file:');
        
        if ~(isnumeric(TFormFilename)&&TFormFilename==0) %check the user has not pressed cancel 
        
            %initialising parameters
            TFormName = [TFormPathname TFormFilename];
            lobject = appData.FilterLobject;
            threshold = appData.localizationThresh;
            
            if (appData.checkboxfirstframeG == 1)
                FirstGreenFrame = 1;
                alternationPeriod = 1;
            else
                FirstGreenFrame = 2;
                alternationPeriod = 2;
            end    
            
            tif = appData.checkboxLoadTifFile;
            alttype = appData.AltType;
            startframe = appData.startframe;
                  
            %loop over loaded movies
            for ii = 1:appData.nFiles
            
            if appData.nFiles == 1
                loadname  = [moviePathname movieFilename];
            else
                loadname = [moviePathname movieFilename{1,ii}];
            end
           
            if appData.checkboxFindFilterParams == 1
                FilterAverParam = findFilterParameters(loadname,appData.checkboxLoadTifFile,appData.findFilterParamsOption);
                appData.localizationThresh = FilterAverParam(1,1);
                appData.FilterLnoise = FilterAverParam(1,2);
                appData.FilterLobject = FilterAverParam(1,3);
            end
   
            %FRETSTORM(loadname,threshold,lobject,FirstGreenFrame,TFormName,alternationPeriod,appData.checkboxLoadTifFile);
            savename = FRETSTORMAPModified(loadname,threshold,lobject,FirstGreenFrame,TFormName,alternationPeriod,tif,alttype,startframe,appData.leftRim);
            end
            
        if appData.nFiles == 1
            newData = importdata(savename);%load position data in workspace
            if alttype == 2
            appData.data = newData(2:3:end,:);%only DA positions (red channel)   
            else
            appData.data = newData(2:2:end,:);%only DA postions (red channel)
            end
            clear newData
        end    
            
        end    
        end
   
    case 'FRETtracking_Callback'    %calls trackWithDummy.m
        
        [dataFilename, dataPathname] =...
            uigetfile('*.out', 'FRET localization data:','MultiSelect', 'on');
        
        if ~(isnumeric(dataFilename)&&dataFilename==0) %check the user has not pressed cancel
        
        info = whos('dataFilename');
        if strcmp(info.class,'char')
            appData.nFiles = 1;
        else
            appData.nFiles = numel(dataFilename);
        end
        
        for ii = 1:appData.nFiles
            
            if appData.nFiles == 1
                loadname  = [dataPathname dataFilename];
                appData.dataPathname = dataPathname;
                appData.dataFilename = dataFilename;
            else
                loadname = [dataPathname dataFilename{1,ii}];
            end
            
            newData = importdata(loadname);
            if isstruct(newData)
                appData.data = newData.data;
            else
                appData.data = newData;
            end
            clear newData;
            
            %sort out pos-files depending on alttype for tracking
            
            colNum = length(appData.data(1,:));
            
            if appData.AltType ==2
                %ALEX
                 % standard indexing red-channel
                 ind = find(mod(1:1:length(appData.data(2:2:end,1)),3)>0);
                 localizationsR = appData.data(ind,:);
                 % standard indexing green-channel
                 localizationsG = appData.data(1:3:end,:);
            else
                %CW
                 % standard indexing red-channel
                 localizationsR = appData.data(2:2:end,:);
                 % standard indexing green-channel
                 localizationsG = appData.data(1:2:end,:);
            end
            %tracking in red-channel
            tracksR = trackWithDummyFRET(localizationsR, appData.trackParams);
            nMoleculesR = max(tracksR(:,4));
            appData.tracks = tracksR;
            appData.nMolecules = nMoleculesR; %#mol in DA-channel
            set(handles.results,'String',['nMolecules DA = ' num2str(nMoleculesR)]);
            %tracking in green-channel
            tracksG = trackWithDummyFRET(localizationsG, appData.trackParams);
            nMoleculesG = max(tracksG(:,4));
                      
            if appData.checkboxSaveTracks == 1
%           choose where to save the file                     
%                 [tracksFilename,tracksPathname]=uiputfile([loadname '.tracks'],'Save tracking file as');
%                 loadnametracks = [tracksPathname tracksFilename];
%                 save([loadnametracks , '_win',num2str(appData.trackParams.maxDisp),'.tracks']);

%           save file in same directory with useful info like tracking window size
            save([loadname(1:end-4),'_win',num2str(appData.trackParams.maxDisp),'_mem',num2str(appData.trackParams.mem),'_.tracks']);    
            end
            
            if appData.checkboxAviMovie == 1
            avimovie(tracksR,appData.checkboxLoadTifFile,appData.plotTracksMinSteps);
            end;
            
            if appData.checkboxTifStack == 1
            tiffstack(tracksR,appData.checkboxLoadTifFile,appData.plotTracksMinSteps);
            end

        end
        
        end
        
    case 'FRETROIs_Callback'    %calls selectFRETROIs, now it calls selectROIs
        alttype=appData.AltType;
        if ~isfield(appData,'data')
            [appData.dataFilename, appData.dataPathname] = uigetfile('*.out', 'MatFile data:','MultiSelect', 'on');
            if ~(isnumeric(appData.dataFilename)&&appData.dataFilename==0) %check the user has not pressed cancel
            newData = importdata([appData.dataPathname appData.dataFilename]);
            if isstruct(newData)
                appData.data = newData.data;
            else
                appData.data = newData;
            end
            clear newData;
            end
        end

        if isfield(appData,'data')
        FRETROIData = selectROIs(appData.data, alttype); %selectROIs(appData.data);
        appData.data = FRETROIData;
        save([appData.dataPathname appData.dataFilename 'FRETROI.out'],'FRETROIData','-ascii')
        end               
     
        
    case 'startframe_Callback'
        
        startframe = str2double(get(handles.startframe,'String'));
        
        if isnan(startframe)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.startframe = startframe;
        
    case 'leftRim_Callback'
        
        leftRim = str2double(get(handles.leftRim,'String'));
        
        if isnan(leftRim)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.leftRim = leftRim;
        
    case 'checkboxfirstframeG_Callback'
        
        appData.checkboxfirstframeG = get(handles.checkboxfirstframeG,'Value');        
        
    case 'AltType_Callback'
        
        appData.AltType = get(handles.AltType,'Value');
        
    case 'plotFRETdata_Callback'    %Calls FRETSTORMplotting
        %LATER - add possibility of plotting localisations and tracks on
        %top of WL image, and then tracking circles Plochowietz style onto
        %FRET movie tiff stack...
        
%         [WLmovieFilename, WLmoviePathname] =...
%             uigetfile({'*.fits';'*.tif'}, 'Select WL data corresponding to analysed file');
%         
%         if ~(isnumeric(WLmovieFilename)&WLmovieFilename==0) %check the user has not pressed cancel
%         
%             WLloadname = [WLmoviePathname WLmovieFilename];
%         
%         
        [FLmovieFilename, FLmoviePathname] =...
            uigetfile({'*.fits';'*.tif'}, 'Select analysed FL movie');
        
        if ~(isnumeric(FLmovieFilename)&&FLmovieFilename==0) %check the user has not pressed cancel
        
            FLloadname = [FLmoviePathname FLmovieFilename];
        
        
        [FRETSTORMFilename, FRETSTORMPathname] =...
            uigetfile({'*.out'}, 'Select FRETSTORM analysed data');
        
        if ~(isnumeric(FRETSTORMFilename)&&FRETSTORMFilename==0) %check the user has not pressed cancel
        
            analysisname = [FRETSTORMPathname FRETSTORMFilename];
            
            
            tif = appData.checkboxLoadTifFile;
            alttype = appData.AltType;
            timeres = appData.dT;
            plotTracksMinSteps = appData.plotTracksMinSteps;
            
            FRETSTORMplotting(FLloadname, FRETSTORMFilename, tif, alttype, timeres, plotTracksMinSteps);
            
        end 
        end        
        
    case 'overlaymovie_Callback' %Calls FREToverlaymovie

        loadnames = uigetfile({'*.fits';'*.tif'}, 'Select fluorescence movies to overlay','MultiSelect', 'on');

        if ~(isnumeric(loadnames)&&loadnames==0) %check the user has not pressed cancel

            TFormfilename = uigetfile('*.mat', 'Select transformation matrix');

        if ~(isnumeric(TFormfilename)&&TFormfilename==0) %check the user has not pressed cancel    

            FirstGreenFrame = appData.checkboxfirstframeG;
            tif = appData.checkboxLoadTifFile;
            alttype = appData.AltType;
            startframe = appData.startframe;

            FREToverlaymovie(loadnames,TFormfilename, FirstGreenFrame, tif, alttype, startframe);

        end    
        end
    
        
    case 'trackingWindow_Callback'
        
        trackingWindow = str2double(get(handles.trackingWindow,'String'));
        
        if isnan(trackingWindow)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        %appData.trackParams.mem = 0;
        appData.trackParams.dim = 2;
        appData.trackParams.good = 0;
        appData.trackParams.quiet = 0;
        appData.trackParams.maxDisp = trackingWindow;
        
        
    case 'FilterLnoise_Callback'
        
        FilterLnoise = str2double(get(handles.FilterLnoise,'String'));
        
        if isnan(FilterLnoise)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.FilterLnoise = FilterLnoise;
        
    case 'FilterLobject_Callback'
        
        FilterLobject = str2double(get(handles.FilterLobject,'String'));
        
        if isnan(FilterLobject)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.FilterLobject = FilterLobject;    
  
        
        
    case 'plotTracks_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
            dataIn = importdata([appData.tracksPathname appData.tracksFilename]);
            appData.tracks = dataIn.appData.tracks;
            end
            clear dataIn;
        end
        
        if  isfield(appData,'tracks')
        set(handles.results,'String',['nMolecules = ' num2str(max(appData.tracks(:,4)))]);
        plotTracksColourCoded(appData.tracks,appData);
        end
        
        
    case 'plotTracksMinSteps_Callback'
        
        plotTracksMinSteps = str2double(get(handles.plotTracksMinSteps,'String'));
        
        if isnan(plotTracksMinSteps)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.plotTracksMinSteps = plotTracksMinSteps;
  
    case 'numberOfSubplots_Callback'
        
        numberOfSubplots = str2double(get(handles.numberOfSubplots,'String'));
        
        if isnan(numberOfSubplots)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.numberOfSubplots = numberOfSubplots;    
        
        
    case 'trackMemory_Callback'
        
        trackMemory = str2double(get(handles.trackMemory,'String'));
        
        if isnan(trackMemory)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.trackParams.mem = trackMemory;
        
    case 'plotTracksColour_Callback'
        
        plotTracksColour = get(handles.plotTracksColour,'String');
        
        appData.plotTracksColour = plotTracksColour;
               
    case 'checkboxSaveTracks_Callback'
        
        appData.checkboxSaveTracks = get(handles.checkboxSaveTracks,'Value');
        
    case 'checkboxAviMovie_Callback'
        
        appData.checkboxAviMovie = get(handles.checkboxAviMovie,'Value');        

    case 'checkboxTifStack_Callback'
        
        appData.checkboxTifStack = get(handles.checkboxTifStack,'Value');
        
    case 'checkboxFindFilterParams_Callback'
        
        appData.checkboxFindFilterParams = get(handles.checkboxFindFilterParams,'Value');
        
    case 'checkboxLoadTifFile_Callback'
        
        appData.checkboxLoadTifFile = get(handles.checkboxLoadTifFile,'Value');
   
    case 'checkboxDallSelect_Callback'
        
        appData.checkboxDallSelect = get(handles.checkboxDallSelect,'Value');    
        
    case 'MSDcurve_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
            appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        
        
        if  isfield(appData,'tracks')
            
        [appData.D appData.locD appData.locFrameD appData.locFrameConfineD] =...
            averageFrameMSD(appData.tracks, appData);
        
        assignin('base', 'MSDresults', [appData.D appData.locD appData.locFrameD appData.locFrameConfineD]);
        
        disp_str = {['D = ' num2str(appData.D)];...
            ['locD = ' num2str(appData.locD)];...
            ['locFrameD = ' num2str(appData.locFrameD)];...
            ['locFrameConfineD = ' num2str(appData.locFrameConfineD)]};
        
        set(handles.results,'String',disp_str);
        
        end
        
        
    case 'pixel_Callback'
        
        pixel = str2double(get(handles.pixel,'String'));
        
        if isnan(pixel)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.pixel = pixel;
        
        
    case 'dT_Callback'
        
        dT = str2double(get(handles.dT,'String'));
        
        if isnan(dT)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.dT = dT;
        
        
    case 'sigmaNoise_Callback'
        
        sigmaNoise = str2double(get(handles.sigmaNoise,'String'))/appData.pixel;
        
        if isnan(sigmaNoise)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.sigmaNoise = sigmaNoise;
        
        
    case 'bacLength_Callback'
        
        bacLength = str2double(get(handles.bacLength,'String'));
        
        if isnan(bacLength)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.bacLength = bacLength;
        appData.bacArea = appData.bacWidth*appData.bacLength / (appData.pixel^2);
        
        
    case 'bacWidth_Callback'
        
        bacWidth = str2double(get(handles.bacWidth,'String'));
        
        if isnan(bacWidth)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.bacWidth = bacWidth;
        %also update bacArea
        appData.bacArea = appData.bacWidth*appData.bacLength / (appData.pixel^2);
        
        
    case 'maxStep_Callback'
        
        maxStep = str2double(get(handles.maxStep,'String'));
        
        if isnan(maxStep)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.maxStep = maxStep;
        
        
    case 'histD_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
            appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        
        if  isfield(appData,'tracks')
        appData.histD = histD(appData.tracks,appData);
        end
        
        
    case 'rangeD_Callback'
        
        rangeD = str2num(get(handles.rangeD,'String'));
        
        if isnan(rangeD)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.rangeD = rangeD;
        
        
    case 'DhistMinSteps_Callback'
        
        DhistMinSteps = str2double(get(handles.DhistMinSteps,'String'));
        
        if isnan(DhistMinSteps)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.DhistMinSteps = DhistMinSteps;
        
        
    case 'CDFcurve_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
            appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        
        if  isfield(appData,'tracks')        
        [appData.CDFresults] =...
            CDFAnalysis(appData.tracks, appData, appData.initGuess, appData.CDFfittingOption);
        
        assignin('base', 'CDFresults', appData.CDFresults);
        
        disp_str = {['CDFresults = ' num2str(appData.CDFresults)]};
        
        set(handles.results,'String',disp_str);
        end
        
        
    case 'CDFfittingOptions_Callback'
        
        appData.CDFfittingOption = get(handles.CDFfittingOptions,'Value');
        
        if appData.CDFfittingOption == 1
            appData.initGuess = [];
            set(handles.initGuess,'String',num2str(appData.initGuess));
            disp_str = 'CDF curve: no fit - no initial guess';
            set(handles.initGuessText,'String',disp_str);
            
        elseif appData.CDFfittingOption == 2
            appData.initGuess = 1;
            set(handles.initGuess,'String',num2str(appData.initGuess));
            disp_str = 'CDF curve: initial guess for D [um^2/s]';
            set(handles.initGuessText,'String',disp_str);
            
        elseif appData.CDFfittingOption == 3
            appData.initGuess = [0.5 0.1 1];
            set(handles.initGuess,'String',...
                [num2str(appData.initGuess(1)) ' ' num2str(appData.initGuess(2)) ' ' num2str(appData.initGuess(3))]);
            disp_str = 'CDF curve: initial guess for a, fixed values for D1, D2 [um^2/s]';
            set(handles.initGuessText,'String',disp_str);
            
        elseif appData.CDFfittingOption == 4
            appData.initGuess = [0.5 0.1 1];
            set(handles.initGuess,'String',...
                [num2str(appData.initGuess(1)) ' ' num2str(appData.initGuess(2)) ' ' num2str(appData.initGuess(3))]);
            disp_str = 'CDF curve: initial guess for a, D1, D2 [um^2/s]';
            set(handles.initGuessText,'String',disp_str);
            
        end

   case 'PlotTracksSelection_Callback'
        
        appData.PlotTracksSelection = get(handles.PlotTracksSelection,'Value');    
        
    case 'TracksSelection_Callback'
        
        appData.TracksSelection = get(handles.TracksSelection,'Value');
    
    case 'findFilterParamsOption_Callback'
        
        appData.findFilterParamsOption = get(handles.findFilterParamsOption,'Value');    
        
    case 'initGuess_Callback'
        
        initGuess = str2num(get(handles.initGuess,'String'));
        
        if isnan(initGuess)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.initGuess = initGuess;
        
        
    case 'plotLocalizations_Callback'
        
        if ~isfield(appData,'data')
            [appData.dataFilename, appData.dataPathname] = uigetfile('*.out', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.dataFilename)&&appData.dataFilename==0) %check the user has not pressed cancel
            newData = importdata([appData.dataPathname appData.dataFilename]);
            if isstruct(newData)
                appData.data = newData.data;
            else
                appData.data = newData;
            end
            clear newData;
            end
        end
        
        if isfield(appData,'data')
            
        % standard indexing
        pos(:,1) = appData.data(:,2);
        pos(:,2) = appData.data(:,3);
        
        if appData.holdFigureCheckbox
            figure(appData.figureHandle);
        else
            figure;
        end
        hold all
        plot(pos(:,1),pos(:,2),'.','MarkerSize',3.0);
        %[X,Y] = meshgrid(pos(:,1),pos(:,2));
        %imagesc(pos(:,4));colorbar;

        xlabel('x [pixels]');
        ylabel('y [pixels]');
        axis image;
        
        end
        
        
    case 'selectROIs_Callback'
        
        if ~isfield(appData,'data')
            [appData.dataFilename, appData.dataPathname] = uigetfile('*.out', 'MatFile data:','MultiSelect', 'on');
            if ~(isnumeric(appData.dataFilename)&&appData.dataFilename==0) %check the user has not pressed cancel
            newData = importdata([appData.dataPathname appData.dataFilename]);
            if isstruct(newData)
                appData.data = newData.data;
            else
                appData.data = newData;
            end
            clear newData;
            end
        end

        if isfield(appData,'data')        
        ROIData = selectROIs(appData.data);
        appData.data.data = ROIData;
        alttype = appData.AltType;
        saveName = [appData.dataPathname appData.dataFilename];
        save([saveName(1:end-4) '_ROI.out'],'ROIData','-ascii');
        end
        
    case 'findBindingDiffusionTracks_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
            appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        
        if  isfield(appData,'tracks')   
        findBindingDiffusionTracks(appData.tracks, appData)
        end
        
        
    case 'DThresh_Callback'
        
        DThresh = str2num(get(handles.DThresh,'String'));
        
        if isnan(DThresh)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.DThresh = DThresh;
        
        
    case 'TwoSpeciesMSDThreshold_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
            appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        
        if  isfield(appData,'tracks')   
        
        [appData.diffusionFraction, appData.D1, appData.D2] =...
            twoSpeciesMSDThreshold(appData.tracks, appData);
        
        assignin('base', 'MSDThreshResults', [appData.diffusionFraction, appData.D1, appData.D2]);
        
        disp_str = {['D histogram threshold results = '...
            num2str(appData.diffusionFraction) ', ' num2str(appData.D1) ', '  num2str(appData.D2)]};
        
        set(handles.results,'String',disp_str);
        
        end
        
        
    case 'DhistThresh_Callback'
        
        DhistThresh = str2num(get(handles.DhistThresh,'String'));
        
        if isnan(DhistThresh)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.DhistThresh = DhistThresh;
        
        
    case 'holdFigureCheckbox_Callback'
        
        appData.holdFigureCheckbox = get(handles.holdFigureCheckbox,'Value');
        if appData.holdFigureCheckbox
            appData.figureHandle = figure;
        end
        
        
    case 'curveColor_Callback'
        
        appData.curveColor = get(handles.curveColor,'String');
        
        
    case 'combineTracks_Callback'
        
        appData.tracks = combineTracks;
        
        
    case 'plotBrightfield_Callback'
        
        %modification AS: leftRim set in FRET mode is subtracted, so FRET
        %localizations can be overlayed with brightfield wo shift
        
        [brightfieldFilename, brightfieldPathname] = ...
            uigetfile({'*.fits';'*.tif'}, 'Fits/Tif data:');
        
        if ~(isnumeric(brightfieldFilename)&&brightfieldFilename==0) %check the user has not pressed cancel        
            if (appData.checkboxLoadTifFile == 0)
        ImageInfo = fits_read_header([brightfieldPathname brightfieldFilename]);
        imageLim = [1 ImageInfo.NAXIS1 1 ImageInfo.NAXIS2];
       
        %modification AS: leftRim set in FRET mode is subtracted, so FRET
        %localizations can be overlayed with brightfield, brightfield
        %flipped for same reason
        imageI = double(flipud(getFrame([brightfieldPathname brightfieldFilename],2)));
            else
            info = imfinfo([brightfieldPathname brightfieldFilename]);
            imageI = imread([brightfieldPathname brightfieldFilename], 2, 'Info', info);
            imageI = flipud(imageI);
            end
        
        imageI=imageI(:,appData.leftRim:end);
        
        if appData.holdFigureCheckbox;
            figure(appData.figureHandle);
            hold all;
        else
            figure;
        end
        
        %imshow(imageI,[min(min(imageI)) max(max(imageI))]);
        imagesc(flipud(imageI)), colormap(gray), axis equal;
        
        end
        
        
    case 'clearData_Callback'
        
        if isfield(appData,'data')
            appData = rmfield(appData,'data');
        end
        if isfield(appData,'tracks')
            appData = rmfield(appData,'tracks');
        end
        set(handles.results,'String','');
        
    case 'UBTimeTrace_Callback'
        
        bindingTime = findBindingDiffusionTracks_moving(appData.tracks,appData,appData.filterDwindow,appData.minImmobilizedFrames);
        
    case 'filterDwindow_Callback'
        
        filterDwindow = str2num(get(handles.filterDwindow,'String'));
        
        if isnan(filterDwindow)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.filterDwindow = filterDwindow;
        
    case 'minImmobilizedFrames_Callback' 
        
        minImmobilizedFrames = str2num(get(handles.minImmobilizedFrames,'String'));
        
        if isnan(minImmobilizedFrames)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.minImmobilizedFrames = minImmobilizedFrames;
        
        
        
    %AS Comovement
    case 'comovement_Callback'
        
        %Load Movie
        [movieFilename, moviePathname] =...
            uigetfile({'*.fits';'*.tif'}, 'Image data:','MultiSelect', 'on');
        
        if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
        
            %load Tform
            [TFormFilename, TFormPathname] =...
                uigetfile('*.mat', 'Transform .mat-file:');
        
            if ~(isnumeric(TFormFilename)&&TFormFilename==0) %check the user has not pressed cancel 
        
                %initialising parameters
                TFormName = [TFormPathname TFormFilename];
                lobject = appData.FilterLobject;
                Redthreshold = appData.redThreshold;
                Greenthreshold = appData.greenThreshold;
            
                if (appData.checkboxfirstframeG == 1)
                    FirstGreenFrame = 1;
                    alternationPeriod = 1;
                else
                    FirstGreenFrame = 2;
                    alternationPeriod = 2;
                end    
            
                tif = appData.checkboxLoadTifFile;
                alttype = appData.AltType;
                startframe = appData.startframe;
                  
                %loop over loaded movies
                for ii = 1:appData.nFiles
            
                    if appData.nFiles == 1
                        loadname  = [moviePathname movieFilename];
                    else
                        loadname = [moviePathname movieFilename{1,ii}];
                    end
           
                    if appData.checkboxFindFilterParams == 1
                        FilterAverParam = findFilterParameters(loadname,appData.checkboxLoadTifFile,appData.findFilterParamsOption);
                        appData.localizationThresh = FilterAverParam(1,1);
                        appData.FilterLnoise = FilterAverParam(1,2);
                        appData.FilterLobject = FilterAverParam(1,3);
                    end
   
                    %FRETSTORM(loadname,threshold,lobject,FirstGreenFrame,TFormName,alternationPeriod,appData.checkboxLoadTifFile);
                    [savenameGreen, savenameRed] = ComovStormAS(loadname,Greenthreshold, Redthreshold,lobject,FirstGreenFrame,TFormName, alternationPeriod, tif,alttype,startframe,appData.leftRim);
                end   
            
            end    
        end
        
        
    case 'comovTrackConnect_Callback' %calls ComovCalc
        %needs to get two tracks one for red one for green, calculate
        %distance according to threshold, output: Co-Moving localizations
        %and tracks
        
        %Load tracks green channel
        
        [greenTrackFilename, greenTrackPathname] = uigetfile('*.tracks', 'Green Tracks:');
        
        if ~(isnumeric(greenTrackFilename)&&greenTrackFilename==0) %check the user has not pressed cancel
        
            %load red tracks
            [redTrackFilename, redTrackPathname] =...
                uigetfile({'*.tracks'}, 'Red Tracks:');
        
            if ~(isnumeric(redTrackFilename)&&redTrackFilename==0) %check the user has not pressed cancel
                
                %load localizations green
                [greenLocalFilename, greenLocalPathname] =...
                    uigetfile({'*.out'}, 'Green Localizations:');
                
                if ~(isnumeric(greenLocalFilename)&&greenLocalFilename==0) %check the user has not pressed cancel
                    % Load red localizations
                    [redLocalFilename, redLocalPathname] =...
                        uigetfile({'*.out'}, 'Red Localizations:');
                    
                    %Load arguments
                    loadGreenLocal  = [greenLocalPathname greenLocalFilename];
                    loadRedLocal  = [redLocalPathname redLocalFilename];
                    loadGreenTrack  = [greenTrackPathname greenTrackFilename];
                    loadRedTrack  = [redTrackPathname redLocalFilename];
                    lobject = appData.FilterLobject;
                    distanceThreshold=4;
                    stepnumber=4;
                    
                    
                    %call ComovCalc
                    [saveTracks, saveLocal ] = ComovCalc( loadGreenLocal, loadRedLocal, loadGreenTrack, loadRedTrack, lobject,distanceThreshold, stepnumber);
                    
                end   
            
            end    
        end
        
        
              
 end

%update handles
setappdata(handles.figure1, 'appData', appData);
guidata(handles.figure1, handles);

