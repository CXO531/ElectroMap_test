function varargout = ElectroMap(varargin)

% Main fucntion for running ElectroMap. 
% Chris O'Shea and Ting Yue Yu, University of Birmingham 
% Version 1.0
% Release Date - 

% ElectroMap MATLAB code for ElectroMap.fig
%      ElectroMap, by itself, creates a new ElectroMap or raises the existing
%      singleton*.
%
%      H = ElectroMap returns the handle to a new ElectroMap or the handle to
%      the existing singleton*.
%
%      ElectroMap('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ElectroMap.M with the given input arguments.
%
%      ElectroMap('Property','Value',...) creates a new ALLM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ElectroMap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      NUstop.  All inputs are passed to ElectroMap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)"
%
% See also: GUIDE, GUIDATA, GUIHANDLESF

% Edit the above text to modify the response to help ElectroMap

% Last Modified by GUIDE v2.5 02-Dec-2018 16:59:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ElectroMap_OpeningFcn, ...
    'gui_OutputFcn',  @ElectroMap_OutputFcn, ...
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


% --- Executes just before ElectroMap is made visible.
function ElectroMap_OpeningFcn(hObject, eventdata, handles, varargin)
% Fucnction for setting values on interface intialisation 
handles.output = hObject;
set(handles.invertopt,'Value',1); %inversion pf signal
set(handles.sfilt,'Value',2); %spatial filtering (gaussian)
set(handles.velout, 'Value', 4); %Velociy outlier removal 

handles.bgon=1; %background on/off switch
handles.bgcol='w'; %backgroung colour

handles.folder_name=[]; %filled when directory folder chosen
handles.fname='opening of GUI hold'; 
handles.lastprocessedfname='opening of GUI hold'; %holds for loaded and processed files

handles.drawcon=0;
handles.conbon=[];
handles.medifilt=1;

masterdir=cd;
if isdeployed == 0
addpath(masterdir,'-frozen');
end

handles.ttpstart=10;
handles.ttpend=90; %time to peak defualt settings 

handles.roinum=1;
handles.roisum=0; %roi defualt settings

set(handles.configure,'Visible','off')

handles.herefromroiload=0; %switch for load roi from .txt file
set(handles.manthresh,'Enable','off') %slider off unless manual trheshold level set
handles.rect=[];
handles.loadedmask=[];

set(handles.resegment,'Enable','off') 
set(handles.B2B,'Enable','off')
handles.herefromsegmentpush=0;
handles.filming = 0; %switch for saving maps to video files

%axis not visible until something is in them
axes(handles.mapaxes); axis off
axes(handles.bgimage);axis off
axes(handles.axes2);axis off
axes(handles.cb); axis off
axes(handles.imageaxes); axis off
zoom xon
handles.isZoomed=0;



handles.fmin=4;handles.fmax=10;handles.fbin=0.05;handles.dfwin=0; %Frequency mapping defaults 

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ElectroMap wait for user response (see UIRESUME)
% uiwait(handles.ElectroMap);


% --- Outputs from this function are returned to the command line.
function varargout = ElectroMap_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushselect.
function pushselect_Callback(hObject, eventdata, handles)
% hObject    handle to pushselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
%% Populate listbox with .mat and .tif files
set(handles.listbox1,'Value',1); %set value to 1 on each new file chosen to stop error
[handles.folder_name]=uigetdir;
workingdir=handles.folder_name;
if isdeployed == 0
cd(workingdir);
end

%get all files in directory
allfiles = dir(handles.folder_name);
[a,b,c]=fileparts(handles.folder_name);
if isdeployed == 0
addpath(handles.folder_name);
end
file_list= {};
count=0;

%Find tif and mat files
for i=1:length(allfiles)
    k = findstr(allfiles(i).name, '.TIF');
    d = findstr(allfiles(i).name, '.tif');
    m = findstr(allfiles(i).name, '.mat');
    l = isdir(allfiles(i).name);
    if isempty(k) ~= 1 && l ~= 1
        count=count+1;
        file=allfiles(i).name;
        file_list{count}=file;
    end
    if isempty(d) ~= 1 && l ~= 1
        count=count+1;
        file=allfiles(i).name;
        file_list{count}=file;
    end
    if isempty(m) ~= 1 && l ~= 1
        count=count+1;
        file=allfiles(i).name;
        file_list{count}=file;
    end
end
set(handles.listbox1,'String',file_list);
guidata(hObject, handles);

% --- Executes on selection change in listbox1.

function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushload.
function pushload_Callback(hObject, eventdata, handles)
% hObject    handle to pushload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
set(handles.listbox2,'Value',1);

%% Get image info from GUI
%image
handles.threshop=get(handles.threshopt,'Value'); %thershold choice
handles.threshman=(get(handles.manthresh,'Value')); %manual thershold setting 
imchoice=get(handles.imagedisp,'Value'); %image to display
cropchoice=get(handles.cropbox,'Value'); %crop image setting 
handles.cropchoice=cropchoice; %0 means no crop, 1 mean new crop, 2 means crop from before

if cropchoice == 1
    handles.rect = [];
end

quinnieopt=get(handles.squareROI,'Value'); %custom ROI setting 

%% file info
chosenfilecontents=cellstr(get(handles.listbox1,'String'));
choice=get(handles.listbox1,'Value');
fname=chosenfilecontents{choice};
if ispc == 1 %change of file setting it mac or pc
handles.fnamenew=[handles.folder_name,'\',fname];
else
handles.fnamenew=[handles.folder_name,'/',fname];
end
tf = strcmp(handles.fname,handles.fnamenew);
if tf == 0
    handles.rect = [];
    handles.images=[];
    handles.lastprocessedfname='Pressing load or changing threshold hold';
end
handles.fname=handles.fnamenew;



%% If custom roi chosen, get rid off manual thresholding
if quinnieopt == 1 || handles.herefromroiload == 1
    handles.threshop = 2;
    handles.threshman = -50000;
    if handles.herefromroiload == 1
        quinnieopt = 0;
    end
end
%% load, crop and reset opt, threshold image
inversion=get(handles.invertopt,'Value');
camopt=0;
axes(handles.imageaxes)
set(handles.resegment,'Enable','off')
set(handles.B2B,'Enable','off')

%% Load new image using OMimload function
if tf == 0
[num_images,handles.newrect,mask,im,handles.I,boundaries,handles.camopt,handles.frame1,handles.fluoim,handles.rois] = OMimload(handles.fname,cropchoice,quinnieopt,handles.threshop,handles.threshman,handles.rect,inversion,camopt,get(handles.imagedisp,'Value'),handles.roinum,handles.roisum);
handles.mask=[];
handles.mask=mask;
handles.im=im; handles.num_images=num_images; handles.rect=handles.newrect;
%handles.loadedmask=[];
%handles.herefromroiload = 0;
if isempty(handles.loadedmask) == 0 && handles.herefromroiload == 1
    handles.mask=[];
    handles.I=[];
    mask=handles.loadedmask;
    mask=uint16(mask);
    boundaries=[];
    boundaries = bwboundaries(mask);
    handles.mask=mask;
    if size(im,1) ~= size(mask,1) || size(im,2) ~= size(mask,2)
        if abs(size(im,1)-size(mask,1)) <= 2 && abs(size(im,2)-size(mask,2)) <= 2
            choice = questdlg('ROI dimensons do not match Image but only slighty off. Would you like to reshape ROI?', ...
                'ROI mismatch', ...
                'Yes','No','Yes');
            switch choice
                case 'Yes'
                [rows,cols]=size(im)
                [rows2,cols2]=size(mask)
                if rows2>rows && cols2>cols
                newmask=zeros(rows,cols);
                for r=1:rows
                    for c=1:cols
                        newmask(r,c)=mask(r,c);
                    end
                end
                end
                if rows>rows2 && cols>cols2
                newmask=zeros(rows,cols)
                for r=1:rows2
                    for c=1:cols2
                        newmask(r,c)=mask(r,c);
                    end
                end
                end
                
                mask=uint16(newmask);
            end
        else
            handles.herefromroiload = 0;
            handles.loadedmask=[];
            guidata(hObject,handles)
            h=errordlg('Loaded ROI dimensions do not match Image')
            waitfor(h)
        end
    end
    handles.I=im.*mask;
    handles.mask=mask;
end
set(handles.cropbox,'Value',0);
handles.boundaries=boundaries;
end

%% Rethreshold loaded image set 
tf
fname
handles.fname
size(handles.images)
pause(2)
if tf == 1
    [num_images,handles.newrect,mask,im,handles.I,boundaries,handles.camopt,handles.frame1,handles.fluoim,handles.rois] = OMimreload(handles.fname,cropchoice,quinnieopt,handles.threshop,handles.threshman,handles.rect,inversion,camopt,get(handles.imagedisp,'Value'),handles.roinum,handles.roisum,handles.frame1,handles.num_images);
    handles.mask=[];
handles.mask=mask;
handles.im=im; handles.num_images=num_images; handles.rect=handles.newrect;
%handles.loadedmask=[];
%handles.herefromroiload = 0;
if isempty(handles.loadedmask) == 0 && handles.herefromroiload == 1
    handles.mask=[];
    handles.I=[];
    mask=handles.loadedmask;
    mask=uint16(mask);
    boundaries=[];
    boundaries = bwboundaries(mask);
    handles.mask=mask;
    if size(im,1) ~= size(mask,1) || size(im,2) ~= size(mask,2)
        if abs(size(im,1)-size(mask,1)) <= 2 && abs(size(im,2)-size(mask,2)) <= 2
            choice = questdlg('ROI dimensons do not match Image but only slighty off. Would you like to reshape ROI?', ...
                'ROI mismatch', ...
                'Yes','No','Yes');
            switch choice
                case 'Yes'
                [rows,cols]=size(im)
                [rows2,cols2]=size(mask)
                if rows2>rows && cols2>cols
                newmask=zeros(rows,cols);
                for r=1:rows
                    for c=1:cols
                        newmask(r,c)=mask(r,c);
                    end
                end
                end
                if rows>rows2 && cols>cols2
                newmask=zeros(rows,cols)
                for r=1:rows2
                    for c=1:cols2
                        newmask(r,c)=mask(r,c);
                    end
                end
                end
                
                mask=uint16(newmask);
            end
        else
            handles.herefromroiload = 0;
            handles.loadedmask=[];
            guidata(hObject,handles)
            h=errordlg('Loaded ROI dimensions do not match Image')
            waitfor(h)
        end
    end
    handles.I=im.*mask;
    handles.mask=mask;
end
set(handles.cropbox,'Value',0);
handles.boundaries=boundaries;
end
%% change pic in GUI
axes(handles.imageaxes);
cla;
if imchoice == 1
    imshow(handles.frame1,[],'InitialMagnification', 400)
    colormap('gray')
    freezeColors
    hold on
    for i=1:size(boundaries,1)
        plot(boundaries{i}(:,2),boundaries{i}(:,1),'r','LineWidth',2);
    end
    hold off
end
if imchoice == 2
    imshow(handles.fluoim,[],'InitialMagnification', 400)
    colormap('jet')
    freezeColors
    hold on
    for i=1:size(boundaries,1)
        plot(boundaries{i}(:,2),boundaries{i}(:,1),'k','LineWidth',2);
    end
    hold off
end



guidata(hObject, handles);

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles = guidata(hObject);
%% get processing options from GUI
%segmentation
minpeakdist = str2num(get(handles.minpeak,'String'));
minpeakdist = ceil(minpeakdist/(1/str2num(get(handles.framerate,'String'))));
handles.minpeakdist=minpeakdist;
segchoice = get(handles.segchoice,'Value');
div=str2num(get(handles.segsize,'String'));
minboundary=str2num(get(handles.minbound,'String'));
minmumofpeaks=str2num(get(handles.minnum,'String'));
handles.avgCL=[];
%Baseline
BLopt=(get(handles.BLopt,'Value'));

%filtering
tfilt=get(handles.tfilt,'Value');
sfilt=get(handles.sfilt,'Value');
sfiltsize=str2num(get(handles.sfiltsize,'String'));

%outlieropts
handles.outlier=get(handles.apdout,'Value');
handles.outliervel=get(handles.velout,'Value');

%inversion
inversion=get(handles.invertopt,'Value');

% frame removal
handles.frameremove=get(handles.removef,'Value');

%% is same file check and process images
chosenfilecontents=cellstr(get(handles.listbox1,'String'))
choice=get(handles.listbox1,'Value');
newfname=chosenfilecontents{choice};
tf = strcmp(newfname,handles.lastprocessedfname);
if tf == 0
    fname=newfname;
    loadnewims=1;
    handles.lastprocessedfname=newfname;
elseif tf == 1 && handles.herefromsegmentpush == 0
    fname=newfname;
    newsettingschoice = questdlg('Re-Process?', ...
        'Re-Load same file', ...
        'Yes','No - Just segment','No - Just segment');
    switch newsettingschoice
        case 'Yes'
            loadnewims = 1;
        case 'No - Just segment'
            loadnewims = 0;
            num_images=handles.num_images;
    end
elseif tf == 1 && handles.herefromsegmentpush == 1
    loadnewims = 0;
end


if loadnewims == 1
    if isempty(handles.rect) == 0
        handles.cropchoice = 1
    end
    handles.averages=[];
    [handles.preimages,images,averages,mask] = OMimprocess(handles.fname,handles.im,handles.rect,handles.num_images,handles.cropchoice,handles.mask,sfilt,sfiltsize,inversion,tfilt,handles.frameremove,handles.camopt,str2num(get(handles.sfiltsigma,'String')));
    handles.waverages=averages;
    handles.images=images;
    handles.averages=averages;
    num_images=handles.num_images;
    handles.mask=mask;
end
set(handles.resegment,'Enable','on')
set(handles.B2B,'Enable','on')
if loadnewims == 0
    averages=handles.averages;
    images=handles.images;
    num_images=handles.num_images;
end


%% Baseline Drift Correction

%Top hat filter 
if BLopt == 1 || BLopt == 4
    th_len=str2num(get(handles.thlen,'String'));
    th_len=(th_len)/str2num(get(handles.framerate,'String'));
    th_len=round(th_len)
    
    se = strel('line', th_len, 0.5);
    BLAV = imopen(averages, se);
    %figure, plot(BL)
end

%Poly 4th degree
if BLopt == 2 || BLopt == 5
    [p,s,mu]=polyfit(1:length(averages),averages,4);
    BLAV=polyval(p,1:length(averages),[],mu);
end
%Poly 11th degree
if BLopt == 3 || BLopt == 6
    [p,s,mu]=polyfit(1:length(averages),averages,11);
    BLAV=polyval(p,1:length(averages),[],mu);
end

% No BL correction 
if BLopt == 7
    BLAV=min(averages)
end

handles.averages = (averages-BLAV); %Baseline subtraction
%% Remove baseline from each pixel

BLAV=BLAV-min(BLAV);
if BLopt == 4 || BLopt == 5 || BLopt == 6
    for t = 1:size(images,3)
        images(:,:,t)=images(:,:,t)+BLAV(t);
    end
    
end

wb=waitbar(0.5,'Removing Baseline');


if BLopt == 1 || BLopt == 2 || BLopt == 3
    for row=1:size(images,1) %%MY BL REMOVAL
        for col=1:size(images,2)
            for frame = 1:size(images,3)
                signal(frame)=images(row,col,frame);
            end
            if inversion == 1
                signal=imcomplement(signal);
            end
            if BLopt == 1
                se = strel('line', th_len, 0.5);
                BL = imopen(signal, se);
            end
            
            if BLopt == 2
                [p,s,mu]=polyfit(1:length(signal),signal,4);
                BL=polyval(p,1:length(images(row,col,:)),[],mu);
            end
            
            if BLopt == 3
                [p,s,mu]=polyfit(1:length(signal),signal,11);
                BL=polyval(p,1:length(signal),[],mu);
            end
            
            if row == 40 && col == 40
               preBL=squeeze(images(40,40,:));
            end
            for frame = 1:size(images,3)
                images(row,col,frame)=images(row,col,frame)+BL(frame);
            end
            images(row,col,:)=images(row,col,:) - min(images(row,col,:)); %make all mins zero
            if row == 40 && col == 40
                for frame = 1:size(images,3)
                    newsignal(frame)=images(row,col,frame);
                end
               preBL=imcomplement(preBL) 
               assignin('base','preBL',preBL-min(preBL))
               assignin('base','BL',(BL-min(BL))')
               newsignal=imcomplement(newsignal)
               newsignal=newsignal-min(newsignal);
               assignin('base','postBL',newsignal')
            end
        end
        
    end
    
end




handles.images=images;
waitbar(0.95,wb,'Segmenting Signal');
wholav=handles.averages;

%% Display Signal
schoice=get(handles.segsignal,'Value');
if schoice == 1
    handles.averages=handles.waverages;
end
if schoice == 2 || schoice == 3 || schoice == 4
if schoice == 2 || schoice == 4
figure,
imshow(handles.frame1, [],'InitialMagnification', 800) 
title('Make your selection and press enter');
[~,rec]=imcrop;
cropfig=gcf;
close(cropfig)
rec=floor(rec);
r1=rec(2)
c1=rec(1)
if r1 == 0 
    r1=1;
end
if c1 == 0 
    c1=1;
end
r2=floor(rec(2)+rec(4));
c2=floor(rec(1)+rec(3));
[ar ac numim]=size(handles.images);
rmask=zeros(ar,ac);
newav=[];
roiim=[];
for r=r1:r2
    for c=c1:c2
         rmask(r,c)=1;   
    end
end
rmask=uint16(rmask);
for j=1:numim
    class(handles.images)
    class(rmask)
    roiim=handles.images(:,:,j).*rmask;
    newav(j)=sum(sum(roiim));
end

% figure,
   newav=imcomplement(newav);
    newav=newav-min(newav);
end
    if schoice == 3 || schoice == 4
        if schoice == 3
            newav=handles.waverages;
        end
  dnewav=smooth(newav)
 dnewav=smooth(diff(dnewav))
dnewav(1:10)=0;
dnewav=dnewav-min(dnewav);
dnewav(1:10)=0;
    newav=[];
    newav=[0,dnewav'];
    end
%save overall average for later
wholav=handles.averages;
handles.averages=newav;
end
set(handles.listbox2,'Value',1)

axes(handles.axes2)
plot(handles.averages)
drawnow()

%% BLremoval
%% Baseline Drift Correction
if schoice == 1 || schoice == 2
    if schoice == 1
        handles.averages=wholav;
    end
if BLopt == 1 || BLopt == 4
    th_len=str2num(get(handles.thlen,'String'));
    th_len=(th_len)/str2num(get(handles.framerate,'String'));
    th_len=round(th_len)
    
    se = strel('line', th_len, 0.5);
    BLAV = imopen(handles.averages, se);
    %figure, plot(BL)
end

if BLopt == 2 || BLopt == 5
    [p,s,mu]=polyfit(1:length(handles.averages),handles.averages,4);
    BLAV=polyval(p,1:length(handles.averages),[],mu);
end

if BLopt == 3 || BLopt == 6
    [p,s,mu]=polyfit(1:length(handles.averages),handles.averages,11);
    BLAV=polyval(p,1:length(handles.averages),[],mu);
end

if BLopt == 7
    BLAV=min(handles.averages)
end

handles.averages = (handles.averages-BLAV);
end
axes(handles.axes2)
plot(handles.averages)
drawnow()

%% DETECT PEAKS
handles.locs=[];handles.q2locs=[];handles.avgCL=[];
[handles.locs,minimas,handles.q2locs,handles.avgCL,handles.numofpeaksoverall,handles.peakheight]=Omseg2...
(handles.averages,str2num(get(handles.peakhigh,'String')),minpeakdist,str2num(get(handles.peakhigh,'String')),minpeakdist,minboundary,segchoice,minmumofpeaks,num_images,div)

%% Zoomed Section
axes(handles.axes2);
origInfo = getappdata(gca, 'matlab_graphics_resetplotview')
handles.isZoomed = 0;
if isempty(origInfo)
    handles.isZoomed = 0;
else
    handles.isZoomed = 1;
end
exposure=1/str2num(get(handles.framerate,'String'));
handles.newlim=get(gca,'XLim')/exposure;
handles.newlim(1)=floor(handles.newlim(1));
handles.newlim(2)=ceil(handles.newlim(2));
%% axes2

handles.avgCL=handles.avgCL.*(1/str2num(get(handles.framerate,'String')));
axes(handles.axes2)
cla
CM =['b','r','g','y','c','m','k'];
exposure=1/str2num(get(handles.framerate,'String'));
handles.averagestime=[0:1:(length(handles.averages)-1)]*exposure;
plot(handles.averagestime, handles.averages,'k'),
xlabel('time (ms) \rightarrow');
ylabel('Fluorescence Intensity');
xlim([0 length(handles.averages)*exposure]);
hold on
handles.locs
minimas
length(handles.averages)
length(handles.averagestime)
plot(handles.averagestime(handles.locs),handles.averages(handles.locs), 'or');
%plot(handles.averagestime(minimas(:,1)),handles.averages(minimas(:,1)), 'xb');
%plot(handles.averagestime(minimas(:,2)),handles.averages(minimas(:,2)), 'ob');
%hy = plot(xlim, [handles.peakhigh handles.peakhigh],'Linestyle', ':', 'Color',[1 0 0]);
before=round(str2num(get(handles.beforeGUI,'String'))/exposure);
after=round(str2num(get(handles.afterGUI,'String'))/exposure);
%if length(handles.locs) > 2
for i = 1:length(handles.q2locs(:,1))
    c=mod(i,6);
    if c == 0;
        c=6;
    end
    handles.q2locs;
    A=(handles.q2locs(i,:));
    if isempty(A) == 1
        errordlg('No constant cycle length regions found. Please adjust pre-process settings')
    end
    if min(A(A>0)) < before %if first peak v.close to beginning it is ignored to stop dim error
        k = find(A);
        A(k(1))=0;
    end
    tstart=min(A(A>0))-before
    if tstart == 0
        tstart = 1;
    end
    tend=max(A)+after;
    if tend > length(handles.averagestime)
        if length(A)>1
        tend = A(end-1)+after;
        elseif length(A) == 1
            tend=length(handles.averagestime);
            c=7;
        end
    end
    tstart
    if tend>length(handles.averagestime);
        tend=length(handles.averagestime)
    end
    tstart
    tend
    plot(handles.averagestime(tstart:tend),handles.averages(tstart:tend),'color',CM(c));
end
%end
colorbar off
hold off

set(gca,'FontSize',8);
ax=gca;
line(get(ax,'XLim'),[handles.peakheight handles.peakheight],'Color','b')
%populate section listbox

section = {};

if handles.numofpeaksoverall == 1
    section{1}=['N/A'];
else
    for i = 1:length(handles.q2locs(:,1))
        length(handles.q2locs(:,1));
        handles.q2locs;
        handles.avgCL;
        section{i}=[num2str(i),' (',num2str((handles.avgCL(2,i))),'ms)'];
    end
end

%% add zoomed section
if handles.isZoomed == 1;
    handles.q2locs
    if length(handles.q2locs(1,:)) < 2 %for single peak seg
        handles.q2locs(:,end+1)=0;
    end
    newline=zeros(1,length(handles.q2locs(1,:)));
    newline(1)=handles.newlim(1);
    newline(2)=handles.newlim(2);
    handles.q2locs=[handles.q2locs;newline];
    newsection='Zoomed Section';
    section{length(section)+1}=newsection;
end

handles.section=section;
%% Couple sections to auto windows
%peak count set to 1 because of first ignored peak, should be changes
handles.q2locs
minimas
peak_count=0;
sec_count=0;
handles.winopt=1;
if length(handles.locs) > 1
if handles.winopt == 2
    for i=1:length(handles.q2locs(:,1))
    sec_count=0;
for j=1:length(handles.q2locs(1,:))
    if handles.q2locs(i,j)~=0
       peak_count=peak_count+1;
       sec_count=sec_count+1;
       i;
       j;
       handles.q2locs;
       minimas;
%       autobefore(sec_count)=handles.q2locs(i,j)-minimas(peak_count,1);
%      autoafter(sec_count)=minimas(peak_count,2)-handles.q2locs(i,j);
autobefore=6;
autoafter=6;
    end
end
autobefore(autobefore>0)
autoafter(autoafter>0)
autobefores(i)=nanmean(autobefore(autobefore>0))*exposure;
autoafters(i)=nanmean(autoafter(autoafter>0))*exposure;
autobefore=[];
autoafter=[];
    end
else
autobefores(i)=20;
autoafters(i)=60;  
end
end


if length(handles.locs) == 1 || length(handles.locs) == 2
if length(handles.locs) == 1
handles.q2locs(1,1)=handles.locs;
autobefores=20;
autoafters=60;
end
if length(handles.locs) == 2
    handles.locs
    handles.q2locs
handles.q2locs=[];
handles.q2locs=handles.locs;
autobefores=20;
autoafters=60;
CL2=handles.locs(2)-handles.locs(1)*exposure;
handles.section
handles.section{1}=[num2str(CL2),'ms']
section=handles.section
end
end
handles.autobefores=round(autobefores)
handles.autoafters=round(autoafters)
beforebuffer=10;
afterbuffer=10;
set(handles.listbox2,'String',section);
axes(handles.mapaxes);
pretty=get(handles.colmap,'String');
jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
delete(wb);

if schoice == 2 || schoice == 3
   handles.averages=wholav;
end

guidata(hObject, handles);

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
axes(handles.axes2)
handles.filming = 0;
exposure=1/str2num(get(handles.framerate,'String'));
handles.averagestime=[0:1:(length(handles.averages(1,:))-1)]*exposure;
plot(handles.averagestime, handles.averages,'k'),
xlabel('time (ms) \rightarrow');
ylabel('Fluorescence Intensity');
xlim([0 length(handles.averages)*exposure]);
hold on
plot(handles.averagestime(handles.locs),handles.averages(handles.locs), 'or');
%hy = plot(xlim, [handles.peakhigh handles.peakhigh],'Linestyle', ':', 'Color',[1 0 0]);
section_choice=get(handles.listbox2,'Value');
beforebuffer=10;
afterbuffer=10;
handles.autobefores
handles.autoafters
% if get(handles.winopt,'Value') == 1
%    set(handles.beforeGUI,'String',num2str(handles.autobefores(section_choice)+beforebuffer));
%    set(handles.afterGUI,'String',num2str(handles.autoafters(section_choice)+afterbuffer));
% end
before=round(str2num(get(handles.beforeGUI,'String'))/exposure);
after=round(str2num(get(handles.afterGUI,'String'))/exposure);

handles.q2locs
section_choice
A=(handles.q2locs(section_choice,:))

if isempty(A) == 1
    errordlg('No constant cycle length regions found. Please adjust pre-process settings')
end
if min(A(A>0)) < before %if first peak v.close to beginning it is ignored to stop dim error
    k = find(A);
    A(k(1))=0;
end
tstart=min(A(A>0))-before;
if tstart == 0
    tstart = 1;
end
tstart
A
tend=max(A)+after
if tend > length(handles.averagestime)
    tend = length(handles.averagestime);
end
plot(handles.averagestime(tstart:tend),handles.averages(tstart:tend),'color','r');
ax=gca;
line(get(ax,'XLim'),[handles.peakheight handles.peakheight],'Color','b')
set(gca,'FontSize',8);


producemaps_Callback(hObject, eventdata, handles)
%guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tfilt.
function segchoice_Callback(hObject, eventdata, handles)
% hObject    handle to tfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tfilt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tfilt


% --- Executes during object creation, after setting all properties.
function segchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function segsize_Callback(hObject, eventdata, handles)
% hObject    handle to segsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of segsize as text
%        str2double(get(hObject,'String')) returns contents of segsize as a double


% --- Executes during object creation, after setting all properties.
function segsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BLopt.
function BLopt_Callback(hObject, eventdata, handles)
% hObject    handle to minpeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minpeak as text
%        str2double(get(hObject,'String')) returns contents of minpeak as a double


% --- Executes during object creation, after setting all properties.
function BLopt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tfilt.
function tfilt_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function tfilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minpeak_Callback(hObject, eventdata, handles)
% hObject    handle to minpeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minpeak as text
%        str2double(get(hObject,'String')) returns contents of minpeak as a double

% --- Executes during object creation, after setting all properties.
function minpeak_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minpeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minnum_Callback(hObject, eventdata, handles)
% hObject    handle to minnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minnum as text
%        str2double(get(hObject,'String')) returns contents of minnum as a double


% --- Executes during object creation, after setting all properties.
function minnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minbound_Callback(hObject, eventdata, handles)
% hObject    handle to minbound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minbound as text
%        str2double(get(hObject,'String')) returns contents of minbound as a double


% --- Executes during object creation, after setting all properties.
function minbound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minbound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in invertopt.
function invertopt_Callback(hObject, eventdata, handles)
% hObject    handle to invertopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invertopt


% --- Executes on selection change in sfilt.
function sfilt_Callback(hObject, eventdata, handles)
% hObject    handle to sfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sfilt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sfilt


% --- Executes during object creation, after setting all properties.
function sfilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sfiltsize_Callback(hObject, eventdata, handles)
% hObject    handle to sfiltsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sfiltsize as text
%        str2double(get(hObject,'String')) returns contents of sfiltsize as a double


% --- Executes during object creation, after setting all properties.
function sfiltsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfiltsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Mapchoice.
function Mapchoice_Callback(hObject, eventdata, handles)
% hObject    handle to Mapchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
contents = cellstr(get(handles.Mapchoice,'String'));
choice=get(handles.Mapchoice,'Value');
axes(handles.bgimage);
if handles.bgon == 0
imshow(handles.frame1,[],'InitialMagnification', 400)
colormap('gray')
freezeColors
else
cla reset
axis off
end
axes(handles.mapaxes);
isochoice=get(handles.isoopt,'Value');
CVmap=[];
if choice == 10
        cla reset
    set(handles.meanchange,'String','');
    set(handles.textchange,'String','');
    t=str2num(get(handles.t,'String'));
   [map,mapSNRdb,alll,allSNRdb]=snrs(handles.averageBeat,handles.mask,10,50,(get(handles.tfilt,'Value')));
    if isempty(map) == 1
        map=0
        alll=0
    end
    dmap=map; %displau map. values altered so displyaed so bg=white etc. Display and expoterd stats NOT based on this map. 
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    %jetcolormap(1,:) = [1,1,1];
    %jetcolormap(end,:) = [0, 0, 0];
    colormap(jetcolormap);
    if get(handles.apdscale,'Value') == 1
         dmap(dmap==0)= NaN;
         him=imshow(dmap,[], 'InitialMagnification', 800);
         set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
         caxis([floor(min(min(alll))) ceil(max(max(alll)))])
    end
    if get(handles.apdscale,'Value') == 2
        dmap(dmap==0)=NaN;
         him=imshow(dmap,[str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))], 'InitialMagnification', 800);
         set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
        caxis([str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))])
    end
    freezeColors
    % colorbar
    axes(handles.cb);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String');
    colormap(colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.apdscale,'Value') == 1
        max(max(alll))
        min(min(alll))
        stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
        hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Signal/Noise';
        axis off
    end
    if get(handles.apdscale,'Value') == 2
        stepp=(str2num(get(handles.cmax,'String'))-str2num(get(handles.cmin,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin,'String')):stepp:str2num(get(handles.cmax,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Signal/Noise';
        axis off
    end
       pos = get(hcb,'Position')
    hcb.Label.Position=[pos(1) pos(2)-1.2];
    stdall=std(alll);
    palll=prctile(alll,[5,50,95]);
    mean(alll)
    stdall
    %pause(10)
    handles.rdata(4,1)=mean(alll);handles.rdata(4,2)=stdall;handles.rdata(4,3)=stdall/sqrt(numel(alll));handles.rdata(4,4)=stdall*stdall;handles.rdata(4,5)=((palll(3)-palll(1))/palll(2));
    
    rownames{1}='APD';rownames{2}='CV';rownames{3}='Amp';rownames{4}='SNR';
    axes(handles.mapaxes);
    set(handles.rtable,'RowName',rownames);
    set(handles.rtable,'data',handles.rdata);
end
if choice == 1 
    cla reset
    set(handles.meanchange,'String','');
    set(handles.textchange,'String','');
    t=str2num(get(handles.t,'String'));
    [map,~,alll]=mapsbaby(get(handles.aptime1,'Value'),str2num(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')),handles.medifilt);
    if isempty(map) == 1
        map=0
        alll=0
    end
    map
    dmap=map; %displau map. values altered so displyaed so bg=white etc. Display and expoterd stats NOT based on this map. 
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    %jetcolormap(1,:) = [1,1,1];
    %jetcolormap(end,:) = [0, 0, 0];
    colormap(jetcolormap);
    if get(handles.apdscale,'Value') == 1
         dmap(dmap==0)= NaN;
         him=imshow(dmap,[], 'InitialMagnification', 800);
         set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
         caxis([floor(min(min(alll))) ceil(max(max(alll)))])
    end
    if get(handles.apdscale,'Value') == 2
        dmap(dmap==0)=NaN;
         him=imshow(dmap,[str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))], 'InitialMagnification', 800);
         set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
        caxis([str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))])
    end
    freezeColors
    % colorbar
    axes(handles.cb);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String');
    colormap(colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.apdscale,'Value') == 1
        max(max(alll))
        min(min(alll))
        stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
        hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Duration (ms)';
        axis off
    end
    if get(handles.apdscale,'Value') == 2
        stepp=(str2num(get(handles.cmax,'String'))-str2num(get(handles.cmin,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin,'String')):stepp:str2num(get(handles.cmax,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Duration (ms)';
        axis off
    end
       pos = get(hcb,'Position')
    hcb.Label.Position=[pos(1) pos(2)-1.2];
    stdall=std(alll);
    palll=prctile(alll,[5,50,95]);
    mean(alll)
    stdall
    %pause(10)
    handles.rdata(1,1)=mean(alll);handles.rdata(1,2)=stdall;handles.rdata(1,3)=stdall/sqrt(numel(alll));handles.rdata(1,4)=stdall*stdall;handles.rdata(1,5)=((palll(3)-palll(1))/palll(2));
    
    [mapSNRr,mapSNRdb,allSNRr,allSNRdb]=snrs(handles.averageBeat,handles.mask,10,50,(get(handles.tfilt,'Value')));
    
     handles.rdata(4,1)=mean(allSNRr);handles.rdata(4,2)=mean(allSNRdb);
%     handles.rdata(4,2)='';handles.rdata(4,3)='';handles.rdata(4,4)='';handles.rdata(4,5)='';
    rownames=get(handles.rtable,'RowName');
    rownames{4}='SNR';
    axes(handles.mapaxes);
    title(['APD', num2str(t)]);
    set(handles.rtable,'RowName',rownames);
    set(handles.rtable,'data',handles.rdata);
end

if choice == 2
    cla
        set(handles.meanchange,'String','');
    set(handles.textchange,'String','');
    if get(handles.actfittimes,'Value') == 1
        [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CVmap] =...
            cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
            0,100,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
    end
    if get(handles.actfittimes,'Value') == 2
        [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CVmap] =...
            cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
            str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
    end
    if isochoice == 1
        mini=0;
        maxi=max(max(map));
    elseif isochoice == 2
        mini=str2num(get(handles.isomin,'String'));
        maxi=str2num(get(handles.isomax,'String'));
    end
    dmap=map;
    dmap(dmap==0)= NaN;
    him=imshow(dmap, [0 maxi], 'InitialMagnification', 800),
    hold on
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    colormap(jetcolormap);
    caxis([mini maxi]);
     set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
    title('Activation Map');
    step=floor(maxi/10);
    if step == 0
        step = 1;
    end
    
    
    %  set(isomap_handle, 'AlphaData', 0.5);
    % add scale to image
    scalelength = 1; % 1mm
    pix=str2num(get(handles.pixelsize,'String'));
    scale = round(scalelength/pix);
    %text(0,1,['CV: ',num2str(CV), 'cm/sec'], 'Units', 'Normalized')
    freezeColors
    
    %colorbar
    axes(handles.cb);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    stepp=(maxi-mini)/5;
    hcb.TickLabels=[mini:stepp:maxi];
    hcb.Ticks=[0.01,0.2:0.2:1];
    hcb.Label.String='Time of activation (ms)';
    axis off
    handles.rdata(4,1)=NaN;handles.rdata(4,2)=NaN;handles.rdata(4,3)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,5)=NaN;
    rownames=get(handles.rtable,'RowName');
    rownames{4}='';
    set(handles.rtable,'RowName',rownames);
    set(handles.meanchange,'String',rownames{4});
    set(handles.rtable,'data',handles.rdata);
end
if choice == 3
    cla
        set(handles.meanchange,'String','');
    set(handles.textchange,'String','');
    if get(handles.actfittimes,'Value') == 1
        [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV] =...
            cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,1,str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
            0,100,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
    end
    if get(handles.actfittimes,'Value') == 2
        [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV] =...
            cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,1,str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
            str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
    end
    if isochoice == 1
        mini=0;
        maxi=max(max(map));
    elseif isochoice == 2
        mini=str2num(get(handles.isomin,'String'));
        maxi=str2num(get(handles.isomax,'String'));
    end
dmap=map;
    dmap(dmap==0)= NaN;
    him=imshow(dmap, [0 maxi], 'InitialMagnification', 800),
    hold on
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    colormap(jetcolormap);
    caxis([mini maxi]);
     set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
    step=floor(maxi/10);
    if step == 0
        step = 1;
    end
    %colorbar('FontSize', 8,'Location','southoutside','YTick',[(mini:step:maxi)]);
    % colorbar('FontSize', 20, 'YTickLabel', {'0ms','5ms', '10ms', '15ms', '20ms'});
    % plot velocity vectors
    % axis xy
    scal = str2num(get(handles.scal,'String'));
    hold on
    quiver(quivers_X,quivers_Y,scal*quivers_vx,scal*quivers_vy,0,'k')
    
    % Compile CV map
    factor=str2num(get(handles.pixelsize,'String'))/10;
    [rows cols] = size(map);
    CVmap=zeros(rows,cols);
CVXmap=zeros(rows,cols);
CVXmap(sub2ind(size(CVXmap),quivers_Y,quivers_X)) = quivers_vx;
sCVXmap=CVXmap.*CVXmap;
CVYmap=zeros(rows,cols);
CVYmap(sub2ind(size(CVYmap),quivers_Y,quivers_X)) = quivers_vy;
sCVYmap=CVYmap.*CVYmap;
%construct cv map
CVmap=sqrt(sCVXmap+sCVYmap);
CVmap=CVmap*factor;

    hold off
    title('Activation Map');
    % saveas(gca, fullfile(directory, ['plako114342_100ms','.png']));
    % hold off
    
    %scales the values to work with out camera frame rate and resolution
    %factor = pix/exposure*100; % converts to cm/sec
    %CV = mean(v)*factor;
    % disp(['mean_cv: ', num2str(CV), 'cm/sec']);
    %
    % text(0,1,['CV: ',num2str(CV), 'cm/sec'], 'Units', 'Normalized')
    freezeColors
    %colorbar
    axes(handles.cb);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    hcb.TicksMode='manual'; hcb.TickLabelsMode='manual';
    stepp=(maxi-mini)/5;
    hcb.TickLabels=[mini:stepp:maxi];
    hcb.Ticks=[0.01,0.2:0.2:1];
    hcb.Label.String='Time of activation (ms)';
    axis off
    handles.rdata(4,1)=NaN;handles.rdata(4,2)=NaN;handles.rdata(4,3)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,5)=NaN;
    rownames=get(handles.rtable,'RowName');
    rownames{4}='';
    set(handles.rtable,'RowName',rownames);
    set(handles.meanchange,'String',rownames{4});
    set(handles.rtable,'data',handles.rdata);
    map=CVmap;
end
if choice == 4
    cla
        set(handles.meanchange,'String','');
    set(handles.textchange,'String','');
    if get(handles.actfittimes,'Value') == 1
        [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,v,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout] =...
            cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
            0,100,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
    end
    if get(handles.actfittimes,'Value') == 2
        [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,v,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout] =...
            cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
            str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
    end
    if isochoice == 1
        mini=0;
        maxi=max(max(map));
    elseif isochoice == 2
        mini=str2num(get(handles.isomin,'String'));
        maxi=str2num(get(handles.isomax,'String'));
    end
    axes(handles.mapaxes)
    dmap=map;
    dmap(dmap==0)= NaN;
    him=imshow(dmap, [0 maxi], 'InitialMagnification', 800),
    hold on
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    colormap(jetcolormap);
    caxis([mini maxi]);
     set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
    step=floor(maxi/10);
    if step == 0
        step = 1;
    end
    %colorbar('FontSize', 8,'Location','southoutside','YTick',[(mini:step:maxi)]);
    % colorbar('FontSize', 20, 'YTickLabel', {'0ms','5ms', '10ms', '15ms', '20ms'});
    % plot velocity vectors
    % axis xy
    scal = str2num(get(handles.scal,'String'));
    hold on
    quiver(quivers_Xout,quivers_Yout,scal*quivers_vxout,scal*quivers_vyout,0,'k')
    hold off
        % Compile CV map
    factor=str2num(get(handles.pixelsize,'String'))/10;
    [rows cols] = size(map);
    CVmap=zeros(rows,cols);
CVXmap=zeros(rows,cols);
CVXmap(sub2ind(size(CVXmap),quivers_Yout,quivers_Xout)) = quivers_vxout;
sCVXmap=CVXmap.*CVXmap;
CVYmap=zeros(rows,cols);
CVYmap(sub2ind(size(CVYmap),quivers_Yout,quivers_Xout)) = quivers_vyout;
sCVYmap=CVYmap.*CVYmap;
%construct cv map
CVmap=sqrt(sCVXmap+sCVYmap);
CVmap=CVmap*factor;
    title('Activation Map');
    % saveas(gca, fullfile(directory, ['plako114342_100ms','.png']));
    % hold off
    
    %scales the values to work with out camera frame rate and resolution
    %factor = pix/exposure*100; % converts to cm/sec
    %CV = mean(v)*factor;
    % disp(['mean_cv: ', num2str(CV), 'cm/sec']);
    %
    % text(0,1,['CV: ',num2str(CV), 'cm/sec'], 'Units', 'Normalized')
    freezeColors
    
    axes(handles.cb);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    stepp=(maxi-mini)/5;
    hcb.TickLabels=[mini:stepp:maxi];
    hcb.Ticks=[0.01,0.2:0.2:1];
    hcb.Label.String='Time of activation (ms)';
    axis off
        handles.rdata(4,1)=NaN;handles.rdata(4,2)=NaN;handles.rdata(4,3)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,5)=NaN;
    rownames=get(handles.rtable,'RowName');
    rownames{4}='';
    set(handles.rtable,'RowName',rownames);
    set(handles.meanchange,'String',rownames{4});
    set(handles.rtable,'data',handles.rdata);
    map=CVmap;
end

if choice == 5
    cla
    axes(handles.mapaxes);
    axis off
    %wb=waitbar(0.9,'Calculating Frequencies');
    tic
    [map]=domfreq(handles.mask,handles.imagerange,str2num(get(handles.framerate,'String')),handles.fmin,handles.fmax,handles.fbin,handles.dfwin);
    toc
    dfs=map(map>0);
    dmap=map;
    if get(handles.apdscale,'Value') == 1
         dmap(dmap==0)= NaN;
         him=imshow(dmap,[min(min(dfs)) max(max(dfs))], 'InitialMagnification', 800);
         set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
         caxis([floor(min(min(dfs))) ceil(max(max(dfs)))])
             pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
          axes(handles.cb);
          
      cla reset
       hcb=colorbar;
             hcb.Location = 'southoutside'
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    hcb.TickLabels=[min(min(dfs)) max(max(dfs))];
    hcb.Ticks=[0 1];
    hcb.Label.String='Dominant Frequency (Hz)';
    end
    
    %delete(wb);
    %imshow(map, [handles.fmin handles.fmax], 'InitialMagnification', 800),
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    %jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    
    if get(handles.apdscale,'Value') == 2
         dmap(dmap==0)= NaN;
         him=imshow(dmap,[str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))], 'InitialMagnification', 800);
         set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
         caxis([str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))])
             pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
          axes(handles.cb);
          
      cla reset
       hcb=colorbar;
             hcb.Location = 'southoutside'
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    hcb.TickLabels=[[str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))]];
    hcb.Ticks=[0 1];
    hcb.Label.String='Dominant Frequency (Hz)';
    end
    freezeColors
    mini=min(min(map(map>0)));
    maxi=max(max(map(map>0)));
    
 
    axis off
    dfs=map(map>0);
    dfsp=prctile(dfs,[5,50,95]);
    handles.rdata(4,1)=mean(dfs);handles.rdata(4,2)=std(dfs);handles.rdata(4,3)=std(dfs)/numel(dfs);handles.rdata(4,4)=std(dfs)*std(dfs);handles.rdata(4,4)=std(dfs)*std(dfs);handles.rdata(4,5)=((dfsp(3)-dfsp(1))/dfsp(2));
    rownames=get(handles.rtable,'RowName');
    rownames{4}='DF';
    set(handles.meanchange,'String',rownames{4});
    set(handles.rtable,'data',handles.rdata);
end
if choice == 6
    wb=waitbar(0.5,'Producing Diastolic Map')
    section_choice=get(handles.listbox2,'Value');
    A=handles.q2locs(section_choice,:);
    A=A(A~=0);
    frame_1=A(1);
    frame_last=A(end);
    exposure = 1/str2num(get(handles.framerate,'String'));
    after=str2num(get(handles.afterGUI,'String'));
    after=round(after/exposure);

    if frame_last+after>handles.num_images
    frame_last=A(end-1)
    end
    [map,~,alll]=DInt(str2num(get(handles.framerate,'String')),handles.I,handles.images,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.afterGUI,'String')),str2num(get(handles.minpeak,'String')),frame_1,frame_last,str2num(get(handles.t,'String')));
    if isempty(alll) == 1
        map=0;
        alll=0;
    end
    numel(alll)
    map(isnan(map)) = 0;
    axes(handles.mapaxes);
    imshow(map,[], 'InitialMagnification', 800);
    title(['Diastolic Interval Distribution']);
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    if get(handles.apdscale,'Value') == 1
        caxis([floor(min(min(alll))) ceil(max(max(alll)))])
    end
    if get(handles.apdscale,'Value') == 2
        caxis([str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))])
    end
    freezeColors
    % colorbar
    axes(handles.cb);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.apdscale,'Value') == 1
        
        stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
        hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Diastolic Interval (ms)';
        axis off
    end
    if get(handles.apdscale,'Value') == 2
        stepp=(str2num(get(handles.cmax,'String'))-str2num(get(handles.cmin,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin,'String')):stepp:str2num(get(handles.cmax,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Diastolic Interval (ms)';
        axis off
        
    end
    dis=alll;
    disper=prctile(dis,[5,50,95]);
    handles.rdata(4,1)=mean(dis);handles.rdata(4,2)=std(dis);handles.rdata(4,3)=std(dis)/numel(dis);handles.rdata(4,4)=std(dis)*std(dis);handles.rdata(4,4)=std(dis)*std(dis);handles.rdata(4,5)=((disper(3)-disper(1))/disper(2));
    rownames=get(handles.rtable,'RowName');
    rownames{4}='DI';
    set(handles.rtable,'RowName',rownames);
    set(handles.rtable,'data',handles.rdata);
    delete(wb)
end
if choice == 7
    cla
    t=str2num(get(handles.t,'String'));
    %[map,~,alll]=ttp2(get(handles.aptime1,'Value'),str2num(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));
    [map,~,alll]=ttpnew(handles.ttpstart,handles.ttpend,str2num(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));
    alll;
    dmap=map;
    title(['Time to Peak Distribution']);
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    colormap(jetcolormap);
    dmap(dmap==0)= NaN;
if get(handles.apdscale,'Value') == 1
         dmap(dmap==0)= NaN;
         him=imshow(dmap,[], 'InitialMagnification', 800);
         set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
         caxis([floor(min(min(alll))) ceil(max(max(alll)))])
    end
    if get(handles.apdscale,'Value') == 2
        dmap(dmap==0)=NaN;
         him=imshow(dmap,[str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))], 'InitialMagnification', 800);
         set(him, 'AlphaData', ~isnan(dmap))
         if handles.bgon == 1
         axis on
         set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', handles.bgcol)
         end
        caxis([str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))])
    end
    freezeColors
    % colorbar
    axes(handles.cb);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.apdscale,'Value') == 1
        max(max(alll))
        min(min(alll))
        stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
        hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Duration (ms)';
        axis off
    end
    if get(handles.apdscale,'Value') == 2
        stepp=(str2num(get(handles.cmax,'String'))-str2num(get(handles.cmin,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin,'String')):stepp:str2num(get(handles.cmax,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Time to Peak (ms)';
        axis off
    end
    ttps=alll;
    ttpsp=prctile(ttps,[5,50,95]);
    handles.rdata(4,1)=mean(ttps);handles.rdata(4,2)=std(ttps);handles.rdata(4,3)=std(ttps)/numel(ttps);handles.rdata(4,4)=std(ttps)*std(ttps);handles.rdata(4,4)=std(ttps)*std(ttps);handles.rdata(4,5)=((ttpsp(3)-ttpsp(1))/ttpsp(2));
    rownames=get(handles.rtable,'RowName');
    rownames{4}='TTP';
    set(handles.rtable,'RowName',rownames);
    set(handles.rtable,'data',handles.rdata);

    
end

if choice == 8
    cla
    wb=waitbar(0.5,'Producing tau map')
    %[map,~,alll,~,~,~,gofmap]=taumap(str2num(get(handles.taustart,'String')),str2num(get(handles.taufinish,'String')),str2num(get(handles.framerate,'String')),handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.r2cut,'String')));
    [map,~,alll,~,~,~,gofmap]=tautest3(str2num(get(handles.taustart,'String')),str2num(get(handles.taufinish,'String')),str2num(get(handles.framerate,'String')),handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.r2cut,'String')));
    map(isnan(map)) = 0;
    cla
    axes(handles.mapaxes)
    imshow(map,[], 'InitialMagnification', 800);
    title('Relaxation Constant Map');
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    if get(handles.apdscale,'Value') == 1
        caxis([floor(min(min(alll))) ceil(max(max(alll)))])
    end
    if get(handles.apdscale,'Value') == 2
        caxis([str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))])
    end
    freezeColors
    % colorbar
    axes(handles.cb);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside';
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    if get(handles.apdscale,'Value') == 1
        max(max(alll));
        min(min(alll));
        stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
        hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Relaxation Constant (ms)';
        axis off
    end
    if get(handles.apdscale,'Value') == 2
        stepp=(str2num(get(handles.cmax,'String'))-str2num(get(handles.cmin,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin,'String')):stepp:str2num(get(handles.cmax,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Relaxation Constant (ms)';
        axis off
    end
    delete(wb)
    taus=alll;
    tausp=prctile(taus,[5,50,95]);
    handles.rdata(4,1)=mean(taus);handles.rdata(4,2)=std(taus);handles.rdata(4,3)=std(taus)/numel(taus);handles.rdata(4,4)=std(taus)*std(taus);handles.rdata(4,4)=std(taus)*std(taus);handles.rdata(4,5)=((tausp(3)-tausp(1))/tausp(2));
    rownames=get(handles.rtable,'RowName');
    rownames{4}='Tau';
    set(handles.rtable,'RowName',rownames);
    set(handles.rtable,'data',handles.rdata);
    
end
if choice == 9
    [map]=fluo_map(str2num(get(handles.framerate,'String')),handles.I,handles.images,get(handles.tfilt,'Value'),handles.averageBeat);
    map(isnan(map)) = 0;
    alll=map(map>0);
    cla
    axes(handles.mapaxes)
    if get(handles.apdscale,'Value') == 1
    imshow(map,[floor(min(min(alll))) ceil(max(max(alll)))], 'InitialMagnification', 800);
    end
    if get(handles.apdscale,'Value') == 2
    imshow(map,[str2num(get(handles.cmin,'String')) str2num(get(handles.cmax,'String'))], 'InitialMagnification', 800);
    end
    title('Amplitude Map');
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);

    freezeColors
    % colorbar
    axes(handles.cb);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside';
    ax = gca;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
        if get(handles.apdscale,'Value') == 1
        stepp=(ceil(max(max(alll)))-floor(min(min(alll))))/5;
        hcb.TickLabels=[floor(min(min(alll))):stepp:ceil(max(max(alll)))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Relaxation Constant (ms)';
        axis off
    end
    if get(handles.apdscale,'Value') == 2
        stepp=(str2num(get(handles.cmax,'String'))-str2num(get(handles.cmin,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin,'String')):stepp:str2num(get(handles.cmax,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Relaxation Constant (ms)';
        axis off
    end
    hcb.Label.String='Signal Level';
    axis off
    handles.rdata(4,1)=NaN;handles.rdata(4,2)=NaN;handles.rdata(4,3)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,4)=NaN;handles.rdata(4,5)=NaN;
    rownames=get(handles.rtable,'RowName');
    rownames{4}='';
    set(handles.meanchange,'String',rownames{4});
    set(handles.rtable,'data',handles.rdata);
end



% Draw Contours
drawcon = handles.drawcon;
if drawcon == 1
    axes(handles.mapaxes);
    hold on
    for j = 0:str2num(handles.conbon):max(max(map))+str2num(handles.conbon);
        mapmask=(map<=j);
        A=map.*mapmask;
        [cons]=bwboundaries(A);
        for i=1:size(cons,1)
            plot(cons{i}(:,2),cons{i}(:,1),'k','LineWidth',1);
        end
    end
end

%    handles.mask=handles.hold_mask; %keep and reinsate overall mask at end
handles.holdmap=map;
handles.holdcvmap=CVmap;
guidata(hObject, handles);
drawnow()

% Hints: contents = cellstr(get(hObject,'String')) returns apdvaluechoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apdvaluechoice


% --- Executes during object creation, after setting all properties.
function Mapchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mapchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mapoptions.
function mapoptions_Callback(hObject, eventdata, handles)
% hObject    handle to mapoptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mapoptions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mapoptions


% --- Executes during object creation, after setting all properties.
function mapoptions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mapoptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ExportMap.
function ExportMap_Callback(hObject, eventdata, handles)
% hObject    handle to ExportMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_fig_children=get(gcf,'children');
Fig_Axes=findobj(GUI_fig_children,'type','Axes');
fig=figure;ax=axes;clf;
new_handle=copyobj(handles.mapaxes,fig);
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])


% --- Executes on button press in actpoints.
function actpoints_Callback(hObject, eventdata, handles)
% hObject    handle to actpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if get(handles.actfittimes,'Value') == 1
    [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,v,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevCV,variCV,SECV] =...
        cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
        0,100,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
end
if get(handles.actfittimes,'Value') == 2
    [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CVv,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevCV,variCV,SECV] =...
        cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
        str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
end
figure, hold on
plot3(act_x,act_y,act_t,'.k')
title('activation points');
zlabel('time (ms)', 'FontSize', 20);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
zlim([0 15])
hold off

% --- Executes on button press in velhist.
function velhist_Callback(hObject, eventdata, handles)
% hObject    handle to velhist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if get(handles.actfittimes,'Value') == 1
    [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,v,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevCV,variCV,SECV] =...
        cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
        0,100,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
end
if get(handles.actfittimes,'Value') == 2
    [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,v,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevCV,variCV,SECV] =...
        cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
        str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
end
figure, histogram(vout,str2num(get(handles.binnumber,'String')));
hold on
%line([mean(vout) mean(vout)],[0 300],'Color','r','Linewidth',3)
xlabel('Conduction Velocity (cm/s)')
ylabel('Number of Pixels')
hold off
% --- Executes on button press in APDdist.
function APDdist_Callback(hObject, eventdata, handles)
% hObject    handle to APDdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
t=str2num(get(handles.t,'String'));
[~,meanapd,alll]=mapsbaby(get(handles.aptime1,'Value'),str2num(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,get(handles.apdout,'Value'),str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));
figure, histogram(alll,str2num(get(handles.binnumber,'String')));
hold on
%line([mean30 mean30],[0 50],'Color','r','Linewidth',3)
xlabel('Action Potential Duration (ms)')
ylabel('Number of Pixels')
hold off

% --- Executes on selection change in apdout.
function apdout_Callback(hObject, eventdata, handles)
% hObject    handle to apdout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns apdout contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apdout


% --- Executes during object creation, after setting all properties.
function apdout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmin_Callback(hObject, eventdata, handles)
% hObject    handle to cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmin as text
%        str2double(get(hObject,'String')) returns contents of cmin as a double


% --- Executes during object creation, after setting all properties.
function cmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmax_Callback(hObject, eventdata, handles)
% hObject    handle to cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmax as text
%        str2double(get(hObject,'String')) returns contents of cmax as a double


% --- Executes during object creation, after setting all properties.
function cmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in velout.
function velout_Callback(hObject, eventdata, handles)
% hObject    handle to velout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.velout, 'Value') == 9
    set(handles.minvel,'String','0.2')
    set(handles.maxvel,'String','0.8')
    set(handles.text64,'String','(0-1)');
elseif get(handles.velout, 'Value') == 2 
    set(handles.minvel,'String','0')
    set(handles.maxvel,'String','100')
    set(handles.text64,'String','cm/s');
end
% Hints: contents = cellstr(get(hObject,'String')) returns velout contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velout


% --- Executes during object creation, after setting all properties.
function velout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minvel_Callback(hObject, eventdata, handles)
% hObject    handle to minvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minvel as text
%        str2double(get(hObject,'String')) returns contents of minvel as a double


% --- Executes during object creation, after setting all properties.
function minvel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxvel_Callback(hObject, eventdata, handles)
% hObject    handle to maxvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxvel as text
%        str2double(get(hObject,'String')) returns contents of maxvel as a double


% --- Executes during object creation, after setting all properties.
function maxvel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushrecal.
function pushrecal_Callback(hObject, eventdata, handles)
% hObject    handle to pushrecal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
contents = cellstr(get(handles.Mapchoice,'String'));
choice=get(handles.Mapchoice,'Value');
axes(handles.mapaxes);
handles.outlier=get(handles.apdout,'Value');
handles.outliervel=get(handles.velout,'Value');
wb=waitbar(0.4,'Calculating APD');
t=str2num(get(handles.t,'String'));
[~,meanapd,~,onedev]=mapsbaby(get(handles.aptime1,'Value'),str2num(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));
waitbar(0.6,wb,'Calculating conduction velocity');
if get(handles.actfittimes,'Value') == 1
    [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,v,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevcv,varicv,SEcv] =...
        cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
        0,100,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
end
if get(handles.actfittimes,'Value') == 2
    [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,v,vout,quivers_Xout,quivers_Yout, quivers_vxout,quivers_vyout,onedevcv,varicv,SEcv] =...
        cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
        str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
end
set(handles.textAPD, 'String', [num2str(meanapd),' +/- ', num2str(onedev), ' ms']);
set(handles.textCV, 'String', [num2str(mean(vout)),' +/- ', num2str(onedevcv), ' cm/s']);




%activation time
tim=act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2num(get(handles.actmax,'String'));
actmin=str2num(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);


if actmax < 100
    Imax=Imax(1);
else Imax=max(tim);
end
Imin=Imin(1);

if actmax < 100
    timmax=Imax*0.01;
else timmax=Imax
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
    timdiff=timmax;
end
if get(handles.checkbox8,'Value')==1
    normfac=225/allpts;
    timdiff=timdiff*normfac;
end

set(handles.actquote, 'String', [num2str(timdiff),' ms']);

delete(wb)
%maps
Mapchoice_Callback(hObject, eventdata, handles)
% hold off
%scales the values to work with out camera frame rate and resolution
%factor = pix/exposure*100; % converts to cm/sec
%CV = mean(v)*factor;
% disp(['mean_cv: ', num2str(CV), 'cm/sec']);
%
% text(0,1,['CV: ',num2str(CV), 'cm/sec'], 'Units', 'Normalized')


guidata(hObject, handles);
% --- Executes on button press in pushmapapply.
function pushmapapply_Callback(hObject, eventdata, handles)
% hObject    handle to pushmapapply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mapchoice_Callback(hObject, eventdata, handles)



% --- Executes on button press in producemaps.
function producemaps_Callback(hObject, eventdata, handles)
% hObject    handle to producemaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
handles.filming=0;
%% Store each peak into array
before=str2num(get(handles.beforeGUI,'String'));
after=str2num(get(handles.afterGUI,'String'));

exposure=1/str2num(get(handles.framerate,'String'));
before = round(before/exposure); %1000 because we are dealing with ms
after = round(after/exposure);
if handles.filming == 0
    wb=waitbar(0.1,'Preparing Images');
    section_choice=get(handles.listbox2,'Value');
end

if handles.filming == 1
    section_choice=handles.filmcount;
end
handles.filming
section_choice
pause(2)
if strcmp(handles.section{section_choice},'Zoomed Section') == 1 && handles.filming ~= 1
    sig=handles.averages(handles.q2locs(section_choice,1):handles.q2locs(section_choice,2));
    %Peak Parameters
    maxfluo = max(handles.averages);
    mo = mode(handles.averages);
    handles.peakheight = maxfluo*str2num(get(handles.peakhigh,'String'));
    minpeakdist = str2num(get(handles.minpeak,'String'));
    minpeakdist = ceil(minpeakdist/(1/str2num(get(handles.framerate,'String'))));
    % find peaks
    [ioioio m] = findpeaks(sig, 'MINPEAKHEIGHT', handles.peakheight, 'MINPEAKDISTANCE', minpeakdist)
    m;
    m=m+handles.q2locs(section_choice,1)-1;
else
    m=handles.q2locs(section_choice,:);
end
f=m(m~=0);
peaks = (length(f)); % ignores last peak as the signal may cut out
% beyond "after" hence matrix dim error
numpeaks = peaks;
timeframe = 1:before+after-1;

% create initial plot
% figure, subplot(1, 2, 1), xlim([0 time(length(timeframe))]);
%figure('name', 'overlaid beats'), subplot(1,2,1);
% xlabel('time (s) \rightarrow');
%ylabel('Fluorescence Intensity');
%hold on
% count = 0;
% oap=[];
% f(1)
% before
% if f(1) <= before
%         startpeak = 2;
% else startpeak = 1;
%     end
% for i = startpeak:(peaks)
%     count = count+1;
%     % stores signals in array
%     start = f(i)-before;
%     finish = f(i)+after;
%     if start < 0
%         start = 1
%     end
%     if length(handles.averages) < finish
%         finish=length(handles.averages);
%     end
%     oap(:,:,i) = handles.averages(start:finish);
% end
%
% % accounts for the for loop starting at index 2
%
% oap = oap(1,:,1:peaks);
%
% % Calculates the signal average (if > 1 peak)
% total = sum(oap,3); % 3 because you want the mean of the 3rd dimension
% averageBeat = total/numpeaks
%
% % calculate the standard deviation of the plot
% stdev = sqrt(sum((std(oap,0,3)).^2));

%% OVERLAYING ALL BEATS TO MAKE AN AVERAGE BEAT
% total action potential duration
APtime = before+after;
[~,~,num]=size(handles.images(:,:,:))
% create empty matrix to fill later on
if APtime <= num
overlay = zeros(size(handles.im,1), size(handles.im,2), APtime);
else
overlay = zeros(size(handles.im,1), size(handles.im,2), num);
end
% skip the first and last AP to forgo any possible errors exceeding matrix
% dimensions
if f(1) <= before
    startloc =2;
else startloc =1;
end
if f(end)+after > num
    endloc=numel(f)-1;
else endloc = numel(f);
end
locRange = startloc:endloc;
if length(handles.q2locs) >1
% for i =1:length(locRange)
%     i
%     f
%     locRange
%     if f(locRange(i)) + after < length(handles.images(1,1,:))
%         locRange2(i)=locRange(i)
%     end
%     locRange=locRange2;
% end

end
% fill matrix
if isempty(locRange)== 1
    errordlg('Peak too close to start/end of file to be analysed with current window settings, next peak analysed')
end

f(locRange(end))+after
if f(locRange(end))+after > length(handles.images(1,1,:))
    locRange=locRange(1:end-1)
end
f
locRange
wsmat=[];
tic

f1=f(locRange);
f1start=(f1(1)-before);
f1end=(f1(end)-after);
handles.imagerange=handles.images(:,:,f1start:f1end);
if length(locRange) > 1 %only need to overlay if more than 1 peak
    wsmat=zeros(size(handles.images,1),size(handles.images,2),numel(-before:after),numel(locRange));
    for x = -before:after
        x+before+1;
        f(locRange)+x;
        numel(f);
        if f(locRange)+after < length(handles.images(1,1,:))
        overlay(:,:,x+before+1) = sum(handles.images(:,:,f(locRange)+x),3)./numel(f);
        overlay(:,:,x+before+1) = overlay(:,:,x+before+1).*double(handles.mask);
        wsmat(:,:,x+before+1,locRange)=handles.images(:,:,f(locRange)+x);
        end
        
    end
end
handles.wsmat=wsmat;
toc

if length(locRange) == 1 %1 beat
    for x = -before:after
        overlay(:,:,x+before+1) = (handles.images(:,:,f(locRange)+x));
        overlay(:,:,x+before+1) = overlay(:,:,x+before+1).*double(handles.mask);
    end
end
handles.cvimages=overlay;
inversion=get(handles.invertopt,'Value');


%% WRITE TO TIFF STACK
% % normalise
%cos = 26/11/16 with new BL removal overlay all negative, so min and max
if handles.numofpeaksoverall > 1
    minI = min(overlay(:))
    maxI = max(overlay(:))
    %
    averageBeat = overlay - minI;
    averageBeat = (2^16-1)*averageBeat./(maxI);
    %
    %make 16 bit
    handles.averageBeat = uint16(averageBeat);
end
if handles.numofpeaksoverall == 1
    disp('hi')
    handles.averageBeat=handles.images;
    %handles.cvimages=handles.images;
end
%% Get numbers
handles.outlier=get(handles.apdout,'Value');
if handles.filming == 0
    waitbar(0.4,wb,'Producing APD map');
    t=str2num(get(handles.t,'String'));
    [apmap,meanapd,~,onedev,var,SE]=mapsbaby(get(handles.aptime1,'Value'),str2num(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')),handles.medifilt);
    
    waitbar(0.6,wb,'Producing Isochronal map');
    
    if get(handles.actfittimes,'Value') == 1
            [actmap,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, onedevCV,varCV,SECV]=...
            cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
            0,200,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
    end
    if get(handles.actfittimes,'Value') == 2
        [actmap,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, onedevCV,varCV,SECV]=...
            cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
            str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
    end
    alll=apmap(apmap>0);
    APp=prctile(alll,[5,50,95]);
    CVp=prctile(vout,[5,50,95]);
    handles.err=[onedev,SE,var,((APp(3)-APp(1))/APp(2))];
    handles.errCV=[onedevCV,SECV,varCV,((CVp(3)-CVp(1))/CVp(2))];

    
    rdata=zeros(4,4);
    rdata(1,1)=meanapd;rdata(1,2)=handles.err(1);rdata(1,3)=handles.err(2);rdata(1,4)=handles.err(3);rdata(1,5)=handles.err(4);
    rdata(2,1)=mean(vout);rdata(2,2)=handles.errCV(1);rdata(2,3)=handles.errCV(2);rdata(2,4)=handles.errCV(3);rdata(2,5)=handles.errCV(4);
    
    [fmap]=fluo_map(str2num(get(handles.framerate,'String')),handles.I,handles.images,get(handles.tfilt,'Value'),handles.averageBeat);
    fmap(isnan(fmap))=0;
    salll=fmap(fmap>0);
    sp=prctile(salll,[5,50,95]);
    rdata(3,1)=mean(salll);
    rdata(3,2)=std(salll);
    rdata(3,3)=std(salll)/sqrt(numel(salll));
    rdata(3,4)=std(salll)*std(salll);
    rdata(3,5)=((sp(3)-sp(1))/sp(2));
    rdata(rdata==0)=NaN;
    handles.rdata=rdata;
    set(handles.rtable,'Data',rdata);
    
    %activation time
    
    tim=actmap(actmap>0);;
    tim=tim-min(tim);
    allpts=numel(tim);
    xbins=0:0.01:max(tim);
    tissueact=100*cumsum(hist(tim,xbins))/allpts
    if isempty(tissueact) ~= 1
    actmax=str2num(get(handles.actmax,'String'));
    actmin=str2num(get(handles.actmin,'String'));
    
    Imax = find(tissueact > actmax);
    Imin = find(tissueact > actmin);
    
    
    if actmax < 100
        Imax=Imax(1);
    else Imax=max(tim);
    end
    Imin=Imin(1);
    
    
    if actmax < 100
        timmax=Imax*0.01;
    else timmax=Imax
    end
    timmin=Imin*0.01;
    timdiff=timmax-timmin;
    if actmin == 0
        timdiff=timmax;
    end
    if get(handles.checkbox8,'Value')==1
        normfac=225/allpts;
        timdiff=timdiff*normfac;
    end
    else timdiff = 0;
    end
    set(handles.actquote, 'String', [num2str(timdiff),' ms']);
    
    d=handles.section
    if strcmp(handles.section{1},'N/A') == 1
        set(handles.CLdisp,'String',['N/A - only one peak']);
    end
    if strcmp(handles.section{1},'N/A') == 0
        if strcmp(handles.section{section_choice},'Zoomed Section') == 1
            set(handles.CLdisp,'String','N/A - Custom Section');
        else
            set(handles.CLdisp,'String',[num2str((handles.avgCL(2,section_choice))),' ms (Frequency = ',num2str(1000/round(handles.avgCL(2,section_choice),-1)),' Hz)']);
        end
    end
    delete(wb)
end
guidata(hObject, handles);
axes(handles.mapaxes);
%% Upadte roi selector


% Make MAPS!!!!!
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
% --- Executes on selection change in threshopt.
function threshopt_Callback(hObject, eventdata, handles)
% hObject    handle to threshopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.threshopt,'Value') == 1
    set(handles.manthresh,'Enable','off')
end
if get(handles.threshopt,'Value') == 2
    set(handles.manthresh,'Enable','on')
end
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns threshopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from threshopt


% --- Executes during object creation, after setting all properties.
function threshopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function manthresh_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
pushload_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function manthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to manthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cropimage.
function cropimage_Callback(hObject, eventdata, handles)
% hObject    handle to cropimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in imagedisp.
function imagedisp_Callback(hObject, eventdata, handles)
% hObject    handle to imagedisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
imchoice=get(handles.imagedisp,'Value');
axes(handles.imageaxes);
cla;
[boundaries] = bwboundaries(handles.mask);

if imchoice == 1
    imshow(handles.frame1,[],'InitialMagnification', 400)
    colormap('gray')
    freezeColors
    hold on
    for i=1:size(boundaries,1)
        plot(boundaries{i}(:,2),boundaries{i}(:,1),'r','LineWidth',2);
    end
    hold off
end

if imchoice == 2
    imshow(handles.fluoim,[],'InitialMagnification', 400)
    colormap('jet')
    freezeColors
       hold on
    for i=1:size(boundaries,1)
        plot(boundaries{i}(:,2),boundaries{i}(:,1),'k','LineWidth',2);
    end
    hold off
end


axes(handles.mapaxes);
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns imagedisp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imagedisp


% --- Executes during object creation, after setting all properties.
function imagedisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imagedisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cropbox.
function cropbox_Callback(hObject, eventdata, handles)
% hObject    handle to cropbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cropbox


% --- Executes on button press in exportvalues.
function exportvalues_Callback(hObject, eventdata, handles)
% hObject    handle to exportvalues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[filename,pathname] = uiputfile({'*.csv';'*.txt';'*.mat'}, 'Save Map values (Isochronal maps with vector overlay can only be saved as .mat files)');
[~,~,ext] = fileparts(filename);
file=[pathname,filename];
choice=get(handles.Mapchoice,'Value');

map=handles.holdmap;
% mat file save
if strcmp('.mat',ext) == 1
    if choice == 1;
        APD_Dist=map;
        save(file,'APD_Dist');
    end
    if choice == 2;
        activation_time=map;
        save(file,'activation_time');
    end
    if choice == 3;
        activation_time=map;
        xpositions=X_pos;
        ypositions=Y_pos;
        xvelocities=X_vel;
        yvelocities=Y_vel;
        velocities=total_vel;
        fractional_up=upmap(handles.cvimages,handles.mask)
        save(file,'activation_time','xpositions','ypositions','xvelocities','yvelocities','velocities','fractional_up');
    end
    if choice == 4;
        activation_time=map;
        xpositions=Xout_pos;
        ypositions=Yout_pos;
        xvelocities=Xout_vel;
        yvelocities=Yout_vel;
        velocities=total_velout;
        fractional_up=upmap(handles.cvimages,handles.mask)
        save(file,'activation_time','xpositions','ypositions','xvelocities','yvelocities','velocities','fractional_up');
    end
    if choice == 6
        DI=handles.holdmap;
        save(file,'DI')
    end
    if choice == 8
        tau=handles.holdmap;
        save(file,'tau')
    end
    if choice == 9
        %RTs
        %[RT30,~,alll]=mapsbaby(get(handles.aptime1,'Value',)(str2num(get(handles.framerate,'String')),30,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));
        [RT50,~,alll]=mapsbaby(get(handles.aptime1,'Value'),str2num(get(handles.framerate,'String')),50,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));
        %[RT80,~,alll]=mapsbaby(get(handles.aptime1,'Value',(str2num(get(handles.framerate,'String')),80,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));
        %RT30over80=RT30./RT80;
        %Tau
        disp('1')
        [~,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, onedevCV,varCV,SECV]=...
            cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
            0,200,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
        cvmean = mean(vout);
        CVdev=onedevCV;
     
        %save(file,'RT30','RT50','RT80','RT30over80','t10BL','t1095','t1090','t1080','t30BL','t3095','t3090','t3080','signal_level','gofmap10BL','gofmap1095','gofmap1090','gofmap1080','gofmap30BL','gofmap3095','gofmap3090','gofmap3080')
        save(file,'RT50','t1080','signal_level','gofmap1080','cvmean','CVdev')
        disp('done')
    end
end

% csv and txt files
if strcmp('.csv',ext) == 1 || strcmp('.txt',ext) == 1
    cho = questdlg('How would you like to export the values?', ...
        'Export', ...
        'Map','List','List');
    switch cho
        case 'Map'
            T=table(map);
            writetable(T,file,'Delimiter',',','WriteVariableNames',false);
            if choice == 2 || choice == 3 || choice == 4
            cvfile=[pathname,'CV_',filename];
            T2=table(handles.holdcvmap);
            writetable(T2,cvfile,'Delimiter',',','WriteVariableNames',false);
            end
        case 'List'
            listmap=reshape(map,numel(map),1);
            listmap=listmap(listmap>0);
            listmap=listmap(isnan(listmap)==0);
            T=table(listmap);
            writetable(T,file,'Delimiter',',','WriteVariableNames',false);
            if choice == 2 || choice == 3 || choice == 4
            cvfile=[pathname,'CV_',filename];
            CVmap=handles.holdcvmap;
            listcvmap=reshape(CVmap,numel(CVmap),1);
            listcvmap=listcvmap(listcvmap>0);
            listcvmap=listcvmap(isnan(listcvmap)==0);
            T=table(listcvmap);
            writetable(T,cvfile,'Delimiter',',','WriteVariableNames',false);
            end
        end

end


% --- Executes on button press in act_movie.
function act_movie_Callback(hObject, eventdata, handles)
% hObject    handle to act_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
meme=0;
gifsave=questdlg('Save File?','Gif Save','Yes','No','Yes');
switch gifsave
    case 'Yes'
        meme=1;
    case 'No'
        meme=0;
end
[map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,~,vout,quivers_Xout,quivers_Yout,quivers_vxout,quivers_vyout, onedevcv]...
    =cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
    str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
if meme == 1
[a,b]=uiputfile('*.gif');
filename=[b,a];
end
h=figure;
hold on
imshow(map,[],'InitialMagnification', 800);
map=map;
maxi=max(max(map)); %ms
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [0, 0, 0];
colormap(jetcolormap);
maxiso=str2num(get(handles.isomax,'String'))
caxis([0 maxiso]);
t3=0;t5=0;t7=0;
delay=0.01;
for j = 1:0.1:ceil(maxi);
    mapmask=(map<j);
    A=map.*mapmask;
    imshow(A, [0 ceil(maxiso)], 'Colormap',jetcolormap, 'InitialMagnification', 400)
    if meme == 1
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if j == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay);
    end
    end
end





function isomin_Callback(hObject, eventdata, handles)
% hObject    handle to isomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of isomin as text
%        str2double(get(hObject,'String')) returns contents of isomin as a double


% --- Executes during object creation, after setting all properties.
function isomin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function isomax_Callback(hObject, eventdata, handles)
% hObject    handle to isomax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of isomax as text
%        str2double(get(hObject,'String')) returns contents of isomax as a double


% --- Executes during object creation, after setting all properties.
function isomax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isomax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in isoopt.
function isoopt_Callback(hObject, eventdata, handles)
% hObject    handle to isoopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns isoopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from isoopt


% --- Executes during object creation, after setting all properties.
function isoopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isoopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in applyiso.
function applyiso_Callback(hObject, eventdata, handles)
% hObject    handle to applyiso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
listbox2_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
axes(handles.mapaxes);
choice = questdlg('Video or Image?', ...
    'Thing', ...
    'Video','Image','Video');
switch choice
    case 'Video'
        [filename,pathname] = uiputfile({'*.avi'}, 'Save .avi image video file of currently displayed maps across all sections' );
        [~,~,ext] = fileparts(filename);
        file=[pathname,filename];
%         if isdepolyed == 0
%         cd(pathname)
%         end
        numsec=length(handles.section);
        vidobj = VideoWriter(file);
        vidobj.FrameRate=1;
        open(vidobj);
        set(gca,'nextplot','replacechildren');
        handles.filming = 1;
        wb=waitbar(0.1,'Producing video file','WindowStyle', 'modal');
        for i=1:numsec
            handles.filmcount = i;
            guidata(hObject, handles);
            Mapchoice_Callback(hObject, eventdata, handles);
            handles = guidata(hObject);
            axes(handles.mapaxes);
            currFrame = getframe;
            writeVideo(vidobj,currFrame);
            waitbar((0.1+0.9*(i/numsec)),wb,'Producing video file');
        end
        close(vidobj);
        delete(wb)
        set(handles.listbox2,'Value',numsec)
        handles.filming = 0;
    case 'Image'
        handles.filming = 1;
        numsec=length(handles.section);
        for i=1:numsec
            handles.filmcount = i;
            handles.b2bimage=1;
            guidata(hObject, handles);
            Mapchoice_Callback(hObject, eventdata, handles);
            GUI_fig_children=get(gcf,'children');
            Fig_Axes=findobj(GUI_fig_children,'type','Axes');
            fig=figure;
            ax=axes;clf;
            new_handle=copyobj(handles.mapaxes,fig);
            pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
            jetcolormap(1,:) = [1, 1, 1];
            colormap(jetcolormap);
            set(gca,'ActivePositionProperty','outerposition')
            set(gca,'Units','normalized')
            set(gca,'OuterPosition',[i/numsec-0.2 i/numsec-0.2 i/numsec i/numsec])
            %set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
            guidata(hObject, handles);
        end
handles.filming==0
end


% --- Executes on button press in segEP.
function segEP_Callback(hObject, eventdata, handles)
% hObject    handle to segEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
segEP
% [filename,pathname] = uiputfile({'*.csv';'*.mat'},'Hi');
% [~,~,ext] = fileparts(filename);
% file=[pathname,filename];
% numsec=length(handles.section);
% wb=waitbar(0,'Producing whole file section analysis: section 1')
% opt = questdlg('Do you want to do a fancy single vector transverse thing??', ...
% 	'Thing', ...
% 	'Yes','No','No');
% switch opt
%     case 'No'
%
% for j=1:numsec
%     handles = guidata(hObject);
% waitbar(j/numsec,wb,['Producing whole file section analysis: section ',num2str(j)]);
% %% Store each peak into array
%
% before=str2num(get(handles.beforeGUI,'String'));
% after=str2num(get(handles.afterGUI,'String'));
% %
% % exposure = time(2); % looks at the time to the first frame this is roughly the exposure
% exposure=1/str2num(get(handles.framerate,'String'));
% before = round(before/exposure) %1000 because we are dealing with ms
%  after = round(after/exposure)
% section_choice = j;
% m=handles.q2locs(section_choice,:) %peak positons of section
% f=m(m~=0);
% peaks = (length(f)); % ignores last peak as the signal may cut out
%                             % beyond "after" hence matrix dim error
% numpeaks = peaks;
% timeframe = 1:before+after+1;
%
% count = 0;
% oap=[];
% if f(1) < before
%         startpeak = 2;
% else startpeak = 1
%     end
% for i = startpeak:(peaks) % 2 to skip the first peak
%     count = count+1;
%     % stores signals in array
%     start = f(i)-before;
%     finish = f(i)+after;
%     oap(:,:,i) = handles.averages(start:finish);
%     %plot(oap(:,:,i),'k'), title(['Overlay of ', num2str(numpeaks),' OAPs']);
% end
%
% % accounts for the for loop starting at index 2
% oap = oap(1,:,1:peaks);
%
% % Calculates the signal average (if > 1 peak)
% total = sum(oap,3); % 3 because you want the mean of the 3rd dimension
% averageBeat = total/numpeaks;
%
% % calculate the standard deviation of the plot
% stdev = sqrt(sum((std(oap,0,3)).^2));
%
% %% OVERLAYING ALL BEATS TO MAKE AN AVERAGE BEAT
% % total action potential duration
% APtime = before+after;
%
% % create empty matrix to fill later on
% overlay = zeros(size(handles.im,1), size(handles.im,2), APtime);
%
% % skip the first and last AP to forgo any possible errors exceeding matrix
% % dimensions
% if f(1) < before
%     startloc =2
% else startloc =1
% end
% locRange = startloc:numel(f);
%
% % fill matrix
% % figure('name', 'overlay of all beats'),
%
% for x = -before:after
%    x+before+1
%    f(locRange)+x
%    numel(f)
%    overlay(:,:,x+before+1) = sum(handles.images(:,:,f(locRange)+x),3)./numel(f);
%    overlay(:,:,x+before+1) = overlay(:,:,x+before+1).*double(handles.mask);
% %    imshow(overlay(:,:,x+before+1),[], 'InitialMagnification', 500);
% %    pause(0.01);
% end
%
% % title('overlay of all beats');
% handles.cvimages=overlay;
% %% WRITE TO TIFF STACK
% % % normalise
% minI = min(overlay(:));
% maxI = max(overlay(:));
% %
% averageBeat = overlay - minI;
% averageBeat = (2^16-1)*averageBeat./(maxI);
% %
% %make 16 bit
% handles.averageBeat = uint16(averageBeat)
% t=str2num(get(handles.t,'String'));
% [~,means(j),~,onedevs(j)]=mapsbaby(get(handles.aptime1,'Value'),str2num(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));
%
%
%     if get(handles.actfittimes,'Value') == 1
%    [~,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, cvonedevs(j)]=...
%     cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
%           0,200,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')));
%     end
%     if get(handles.actfittimes,'Value') == 2
%     [~,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, cvonedevs(j)]=...
%     cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
%           str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')));
%     end
%
%
%
% vouts(j)=mean(vout);
% tim=act_t;
% tim=tim-min(tim);
% allpts=numel(tim);
% xbins=0:0.01:max(tim);
% tissueact=100*cumsum(hist(tim,xbins))/allpts;
%
% actmax=str2num(get(handles.actmax,'String'));
% actmin=str2num(get(handles.actmin,'String'));
%
% Imax = find(tissueact > actmax);
% Imin = find(tissueact > actmin);
%
% if isempty(Imax) == 1 || isempty(Imin) == 1
% Imax=10000
% Imin = 10000
% end
% if actmax < 100
% Imax=Imax(1);
% else Imax=max(tim);
% end
% Imin=Imin(1);
%
% if isempty(Imax) == 1 || isempty(Imin) == 1
% Imax=10000
% Imin = 10000
% end
%
% timmax(j)=Imax*0.01;
% timmin(j)=Imin*0.01;
% timdiff=timmax-timmin;
% if actmin == 0
%     timdiff=timmax;
% end
% if get(handles.checkbox8,'Value')==1
%     normfac=225/allpts;
%     timdiff=timdiff*normfac;
% end
%
% end
% Section=[1:numsec]';
% Cycle_Length=[handles.avgCL(2,1:numsec)]';
% APD=means';
% st_dev=onedevs';
% Conduction_Velocity=vouts';
% st_dev2=cvonedevs';
% Activation_max=timmax';
% Activation_min=timmin';
%
% if strcmp('.csv',ext) == 1
% T=table(Section,Cycle_Length,APD,st_dev,Conduction_Velocity,st_dev2,Activation_max,Activation_min);
% writetable(T,file,'Delimiter',',','WriteVariableNames',true);
% end
% if strcmp('.mat',ext) == 1
%     save(file,'Section','Cycle Length','APD','st_dev','Conduction_Velocity','cv_st_dev');
% end
% delete(wb)
%     case 'Yes'
%          wholetable=zeros(721,numsec);
%          wholetable(:,1)=(-360:1:360);
%         for j=1:numsec
%     handles = guidata(hObject);
% waitbar(j/numsec,wb,['Producing whole file section analysis: section ',num2str(j)]);
% %% Store each peak into array
%
% before=str2num(get(handles.beforeGUI,'String'));
% after=str2num(get(handles.afterGUI,'String'));
% %
% % exposure = time(2); % looks at the time to the first frame this is roughly the exposure
% exposure=1/str2num(get(handles.framerate,'String'));
% before = round(before/exposure) %1000 because we are dealing with ms
%  after = round(after/exposure)
% section_choice = j;
% m=handles.q2locs(section_choice,:) %peak positons of section
% f=m(m~=0);
% peaks = (length(f)); % ignores last peak as the signal may cut out
%                             % beyond "after" hence matrix dim error
% numpeaks = peaks;
% timeframe = 1:before+after+1;
%
% count = 0;
% oap=[];
% if f(1) < before
%         startpeak = 2;
% else startpeak = 1
%     end
% for i = startpeak:(peaks) % 2 to skip the first peak
%     count = count+1;
%     % stores signals in array
%     start = f(i)-before;
%     finish = f(i)+after;
%     oap(:,:,i) = handles.averages(start:finish);
%     %plot(oap(:,:,i),'k'), title(['Overlay of ', num2str(numpeaks),' OAPs']);
% end
%
% % accounts for the for loop starting at index 2
% oap = oap(1,:,1:peaks);
%
% % Calculates the signal average (if > 1 peak)
% total = sum(oap,3); % 3 because you want the mean of the 3rd dimension
% averageBeat = total/numpeaks;
%
% % calculate the standard deviation of the plot
% stdev = sqrt(sum((std(oap,0,3)).^2));
%
% %% OVERLAYING ALL BEATS TO MAKE AN AVERAGE BEAT
% % total action potential duration
% APtime = before+after;
%
% % create empty matrix to fill later on
% overlay = zeros(size(handles.im,1), size(handles.im,2), APtime);
%
% % skip the first and last AP to forgo any possible errors exceeding matrix
% % dimensions
% if f(1) < before
%     startloc =2
% else startloc =1
% end
% locRange = startloc:numel(f);
%
% % fill matrix
% % figure('name', 'overlay of all beats'),
%
% for x = -before:after
%    x+before+1
%    f(locRange)+x
%    numel(f)
%    overlay(:,:,x+before+1) = sum(handles.images(:,:,f(locRange)+x),3)./numel(f);
%    overlay(:,:,x+before+1) = overlay(:,:,x+before+1).*double(handles.mask);
% %    imshow(overlay(:,:,x+before+1),[], 'InitialMagnification', 500);
% %    pause(0.01);
% end
%
%
%
% % title('overlay of all beats');
% handles.cvimages=overlay;
% %% WRITE TO TIFF STACK
% % % normalise
% minI = min(overlay(:));
% maxI = max(overlay(:));
% %
% averageBeat = overlay - minI;
% averageBeat = (2^16-1)*averageBeat./(maxI);
% %
% %make 16 bit
% handles.averageBeat = uint16(averageBeat)
% t=str2num(get(handles.t,'String'));
% [~,means(j),~,onedevs(j)]=mapsbaby(get(handles.aptime1,'Value')(str2num(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));
%
%
%     if get(handles.actfittimes,'Value') == 1
%    [map,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, cvonedevs(j)]=...
%     cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
%           0,200,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')));
%     end
%     if get(handles.actfittimes,'Value') == 2
%     [map,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, cvonedevs(j)]=...
%     cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
%           str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')));
%     end
%     veldistance=.5;%???????????????????????
% veldistancepix=round(veldistance*10000/(str2num(get(handles.pixelsize,'String'))));
%
% m=min(map(map>0));
% [rows cols] = size(map);
% count=0;
% for r=1:rows
%     for c=1:cols
%         if map(r,c) == m
%             count=count+1;
%             minind(count,1)=r;
%             minind(count,2)=c;
%         end
%     end
% end
% rrow=round(mean(minind(:,1)));
% rcol=round(mean(minind(:,2)));
%
% th = 0:pi/180:2*pi;
% xunit = veldistancepix * cos(th) + rcol;
% yunit = veldistancepix * sin(th) + rrow;
% xco=100000000000;
% yco=100000000000;
% count=0;
% mat=[];
% for i =1:numel(xunit)
%
%     time_A=(map(rrow,rcol));
%     if  round(yunit(i)) <= rows &&  round(xunit(i)) <= cols && round(yunit(i)) > 0 &&  round(xunit(i)) > 0
%     time_B=(map(round(yunit(i)),round(xunit(i)))); %y-x switch due to image vs axis difference in matlab
%     else time_B =(time_A -1)
%     end
%     p1=[rcol,rrow];
%     p2=[round(xunit(i)),round(yunit(i))];
%     tim=time_B-time_A;
%     if tim > 0 && (p1(1)-p2(1)) ~= 0 && (p1(2)-p2(2)) ~=0;
%         xcn=p1(1)-p2(1)
%         ycn=p1(2)-p2(2)
%         %if xcn ~= xco
%         %    if ycn ~= yco %if statment to stop repeating points with slightly diffrenct angles from th(i)
%         dc=sqrt((xcn*xcn)+(ycn*ycn))
%         count=count+1;
%         mat(count,1)=p1(1)-p2(1);
%         mat(count,2)=p1(2)-p2(2);
%         mat(count,3)=tim; %time diffrence
%         mat(count,4)=(dc*0.0001*(str2num(get(handles.pixelsize,'String'))))/tim*1000; %speed in cm/s
%         mat(count,5)=th(i)*(180/pi);
%         mat(count,6)=dc; %pixel distance
%         xco=xcn;
%         yco=ycn;
%          %   end
%         %end
%     end
% end
%
% mat=unique(mat,'rows');%get rid of repeat points
% [vel_slow, ind_slow] = min(mat(:,4));
%
% zeroangle=[mat(ind_slow,5)];
% vangle(:,2)=mat(:,4);
% vangle(:,1)=(mat(:,5)-zeroangle);
% mat
%
% vangle = sortrows(vangle,1);
% for i=1:length(vangle(:,1))-1
%     if(vangle(i+1,1)) == vangle (i,1)
%         vangle(i,1)=375 %get rid of repeat angles
%     end
% end
% dcount=1
% angind=[];
% for d=1:1:721
% angind=find(vangle(:,1) == (d-361))
% if isempty(angind) == 0;
%     wholetable(d,j+1)=vangle(angind,2);
%     dcount=dcount+1;
% else wholetable(d,j+1)=NaN;
% end
% angind=[];
% end
% vangle
% vangle=[];
%         end
%         wholetable
%         delete(wb)
%         T=table(wholetable);
%         if strcmp('.csv',ext) == 1
% writetable(T,file,'Delimiter',',','WriteVariableNames',false);
%         end
%
% end


% --- Executes on button press in fold.
function fold_Callback(hObject, eventdata, handles)
% hObject    handle to fold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function framerate_Callback(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framerate as text
%        str2double(get(hObject,'String')) returns contents of framerate as a double


% --- Executes during object creation, after setting all properties.
function framerate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pixelsize_Callback(hObject, eventdata, handles)
% hObject    handle to pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixelsize as text
%        str2double(get(hObject,'String')) returns contents of pixelsize as a double


% --- Executes during object creation, after setting all properties.
function pixelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function binnumber_Callback(hObject, eventdata, handles)
% hObject    handle to binnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binnumber as text
%        str2double(get(hObject,'String')) returns contents of binnumber as a double


% --- Executes during object creation, after setting all properties.
function binnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getpixelinfo.
function getpixelinfo_Callback(hObject, eventdata, handles)
% hObject    handle to getpixelinfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

compare
guidata(hObject, handles);
% --- Executes on button press in phasemap.
function phasemap_Callback(hObject, eventdata, handles)
% hObject    handle to phasemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
phasemapping
%% get refrence signal
[~,~,num]=size(handles.images(:,:,:))
GUI_fig_children=get(gcf,'children');
Fig_Axes=findobj(GUI_fig_children,'type','Axes');
fig=figure,
ax=axes;clf;
new_handle=copyobj(handles.mapaxes,fig);
pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
jetcolormap(1,:) = [1, 1, 1];
colormap(jetcolormap);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
%set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
[rcol rrow] = getpts(gcf);
loc = int32([rrow rcol]);
if size(loc)>1
    loc = [loc(1,1) loc(1,2)];
end
rrow=floor(rrow);
rcol=floor(rcol);
for i=1:num
    refpixelsignal(i)=handles.images(rrow,rcol,i);
end
atimage=gcf;
tfilt=get(handles.tfilt,'Value');
if tfilt == 2
    refpixelsignal = sgolayfilt(refpixelsignal, 3,11);
end
refsig=refpixelsignal;
refsig=imcomplement(refsig);
refsig=refsig-min(refsig);
maxfluo=max(handles.averages);
peakheight = maxfluo*str2num(get(handles.peakhigh,'String'));
[pks locs] = findpeaks(refsig, 'MINPEAKHEIGHT', peakheight, 'MINPEAKDISTANCE', 50);
figure, plot((refsig));
hold on
plot(locs,pks,'or');
refplot=gcf;
choice = questdlg('Are you happy with this refrence signal?', ...
    'Refrence Signal', ...
    'Yes','No','Yes');
switch choice
    case 'Yes'
        close(refplot);
        close(atimage);
        %% do phase
        
        CL=81;
        B = 1/CL*ones(CL,1);
        out = filter(B,1,refpixelsignal);
        shiftedsignal=refpixelsignal-out;
        shiftedsignal=shiftedsignal(CL:length(shiftedsignal));
        %figure, plot(shiftedsignal); title('Shiftedsignal');
        %mean(shiftedsignal)
        hsig=hilbert(shiftedsignal)
        figure, hold on
        axes(gcf)
        x=real(hsig);
        y=imag(hsig);
        a=angle(hsig);
        subplot(2,1,1)
        plot(x,y,'rx')
        subplot(2,1,2)
        plot(1:length(shiftedsignal),a)
        yticks([-pi 0 pi ])
        yticklabels({'-\pi','0','\pi'})
        phasemap(handles.images, handles.mask,refpixelsignal,str2num(get(handles.framerate,'String')),str2num(get(handles.minpeak,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.afterGUI,'String')))
        
    case 'No'
        close(refplot);
        close(atimage);
        h = errordlg('Please reselect refrence point by re-pressing phase map');
end




% --- Executes on button press in compare.
function compare_Callback(hObject, eventdata, handles)
% hObject    handle to compare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
% if handles.comapreinstruct == 0
%     waitfor(msgbox('Select up to 6 pixels (single click) and press Enter'));
%     handles.comapreinstruct = 1;
% end
% GUI_fig_children=get(gcf,'children');
% Fig_Axes=findobj(GUI_fig_children,'type','Axes');
% fig=figure,
% ax=axes;clf;
% new_handle=copyobj(handles.mapaxes,fig);
% jetcolormap = (colormap('jet'));
% jetcolormap(1,:) = [1, 1, 1];
% colormap(jetcolormap);
% set(gca,'ActivePositionProperty','outerposition')
% set(gca,'Units','normalized')
% set(gca,'OuterPosition',[0 0 1 1])
% [rcol rrow] = getpts(gcf);
% loc = int32([rrow rcol]);
% if size(loc)>1
%     loc = [loc(1,1) loc(1,2)];
% end
% rrow=floor(rrow);
% rcol=floor(rcol);
% cm=['b','r','g','m','c','k'];
% [~,~,num] = size(handles.images(:,:,:));
% pixelsignal=zeros(num,length(rrow));
% %Signals
% for j=1:length(rrow)
% for i=1:num
%     pixelsignal(i,j)=handles.images(rrow(j),rcol(j),i);
% end
% inversion=get(handles.invertopt, 'Value');
% if inversion == 1
%     pixelsignal(:,j)=imcomplement(pixelsignal(:,j));
% end
% mini=min(pixelsignal(:,j));
% zeropixelsignal(:,j)=pixelsignal(:,j)-mini;
% tfilt=get(handles.tfilt, 'Value');
% if tfilt == 2
%     zeropixelsignal(:,j)=sgolayfilt(zeropixelsignal(:,j), 3,11);
% end
% if tfilt == 3
%    d=designfilt('lowpassiir', 'PassbandFrequency', 100,'StopbandFrequency', 350, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);
%     zeropixelsignal(:,j)=filtfilt(d,zeropixelsignal(:,j));
% end
% end
% close(fig);
% figure('units','normalized','outerposition',[0 0 1 1])
%
% %image
% subplot(2,3,1)
% imshow(handles.I, [], 'InitialMagnification', 800),
%
% for j=1:length(rrow)
% %positions
% subplot(2,3,1)
% hold on
% plot(rcol(j),rrow(j),'+','MarkerSize',20,'LineWidth',5,'color',cm(j));
% hold off
% %signals
% subplot(2,3,[2:3])
% hold on
% plot(zeropixelsignal(:,j),'LineWidth',2,'color',cm(j))
% xlabel('Time (ms)')
% ylabel('Fluorescence Intensity (arb units)');
% hold off
% axis tight
%
% %APDS
% subplot(2,3,5)
% hold on
% mapsbaby(get(handles.aptime1,'Value')onepix(str2num(get(handles.framerate,'String')),30,handles.mask,handles.images,handles.averageBeat,rrow(j),rcol(j),cm(j),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.afterGUI,'String')));
% xlabel('Time (ms)')
% ylabel('Fluorescence Intensity (arb units)');
% axis tight
% hold off
% end
% guidata(hObject, handles);

pixelinfo

% --- Executes on selection change in errorchoice.
function errorchoice_Callback(hObject, eventdata, handles)
% hObject    handle to errorchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
t=str2num(get(handles.t,'String'));
[~,meanapd,~,onedev,var,SE]=mapsbaby(get(handles.aptime1,'Value'),str2num(get(handles.framerate,'String')),t,handles.I,handles.images,handles.averageBeat,handles.outlier,str2num(get(handles.cmin,'String')),str2num(get(handles.cmax,'String')),get(handles.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')));

if get(handles.actfittimes,'Value') == 1
    [~,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~, onedevCV, varCV,SECV]=...
        cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
        0,200,str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
end
if get(handles.actfittimes,'Value') == 2
    [~,~,~,act_t,~,~,~,~,~,~,vout,~,~,~,~,onedevCV, varCV, SECV]=...
        cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
        str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
end

handles.errCV=[onedevCV,varCV,SECV];
set(handles.textAPD, 'String', [num2str(meanapd),' +/- ', num2str(handles.err(get(handles.errorchoice,'Value'))), ' ms']);

set(handles.textCV, 'String', [num2str(mean(vout)),' +/- ', num2str(handles.errCV(get(handles.errorchoice,'Value'))), ' cm/s']);

tim=act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2num(get(handles.actmax,'String'));
actmin=str2num(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);



if actmax < 100
    Imax=Imax(1);
else Imax=max(tim);
end
Imin=Imin(1);


if actmax < 100
    timmax=Imax*0.01;
else timmax=Imax
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
    timdiff=timmax;
end
if get(handles.checkbox8,'Value')==1
    normfac=225/allpts;
    timdiff=timdiff*normfac;
end
set(handles.actquote, 'String', [num2str(timdiff),' ms']);



guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns errorchoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from errorchoice


% --- Executes during object creation, after setting all properties.
function errorchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to errorchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
conduction


% choice = questdlg('Activation Time Choice', ...
% 	'Refrence Signal', ...
% 	'Free Choice','Set Distance','Activation Curve','Free Choice');
% switch choice
%     case 'Free Choice'
%         close(refplot);
%         close(atimage);
% [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,~,vout,quivers_Xout,quivers_Yout,quivers_vxout,quivers_vyout, onedevcv]...
% = cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'));
% m=min(map(map>0));
% [rows cols] = size(map);
% count=0;
% for r=1:rows
%     for c=1:cols
%         if map(r,c) == m
%             count=count+1;
%             minind(count,1)=r;
%             minind(count,2)=c;
%         end
%     end
% end
% actrow=round(mean(minind(:,1)));
% actcol=round(mean(minind(:,2)));
%
%
% acttime=figure,
% isochoice=get(handles.isoopt,'Value');
% if isochoice == 1
% mini=0;
% maxi=max(max(map));
% elseif isochoice == 2
% mini=str2num(get(handles.isomin,'String'));
% maxi=str2num(get(handles.isomax,'String'));
% end
%
% imshow(map, [0 maxi], 'InitialMagnification', 800),
% hold on
% jetcolormap = (colormap('jet'));
% jetcolormap(1,:) = [1, 1, 1];
% colormap(jetcolormap);
% caxis([mini maxi]);
% freezeColors
% hold on
% plot(actcol,actrow,'r+','MarkerSize',20,'LineWidth',3);
% hold off
%
% [rcol rrow] = getpts(gcf)
% rrow=round(rrow)
% rcol=round(rcol)
% if rrow(1) == 0
%     rrow(1) = 1
% end
% if rrow(2) == 0
%     rrow(2) = 1
% end
% if rcol(1) == 0
%     rcol(1) = 1
% end
% if rcol(2) == 0
%     rcol(2) = 1
% end
% time_A=(map(rrow(1),rcol(1)));
% time_B=(map(rrow(2),rcol(2))); %times in ms
% dist=sqrt(((rrow(2)-rrow(1))^2)+((rcol(2)-rcol(1))^2))*str2num(get(handles.pixelsize,'String')); %distance in um
% speed=dist*0.0001/(time_B-time_A)*1000;
% speed=sqrt(speed*speed);
% time=sqrt((time_B-time_A)^2);
% close(acttime)
% figure,
% imshow(map, [0 maxi], 'InitialMagnification', 800),
% hold on
% jetcolormap = (colormap('jet'));
% jetcolormap(1,:) = [1, 1, 1];
% colormap(jetcolormap);
% caxis([mini maxi]);
% freezeColors
% hold on
% plot(actcol,actrow,'r+','MarkerSize',20,'LineWidth',3);
% p1=[rcol(1),rrow(1)];
% p2=[rcol(2),rrow(2)];
% dp=p2-p1;
% q=quiver(p1(1),p1(2),dp(1),dp(2),0);
% q.Color='k';
% q.LineWidth=3;
% q.MaxHeadSize=0.5;
% text(round(cols*0.8),3,['Speed = ',num2str(speed),' cm/s']);
% text(round(cols*0.8),8, ['Distance = ',num2str(dist/10000),' cm']);
% text(round(cols*0.8),13, ['Time = ', num2str(time), ' ms']);
% hold off
%
%  case 'Set Distance'
%      prompt = {'Enter Distance (cm)'};
% dlg_title = 'Distance Input';
% num_lines = 1;
% defaultans = {'0.1',};
% answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
% veldistance=str2num(answer{1});
% veldistancepix=round(veldistance*10000/(str2num(get(handles.pixelsize,'String'))));
%
% [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,~,vout,quivers_Xout,quivers_Yout,quivers_vxout,quivers_vyout, onedevcv]...
% = cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'));
% m=min(map(map>0));
% [rows cols] = size(map);
% count=0;
% for r=1:rows
%     for c=1:cols
%         if map(r,c) == m
%             count=count+1;
%             minind(count,1)=r;
%             minind(count,2)=c;
%         end
%     end
% end
% actrow=round(mean(minind(:,1)));
% actcol=round(mean(minind(:,2)));
%
%
% acttime=figure,
% isochoice=get(handles.isoopt,'Value');
% if isochoice == 1
% mini=0;
% maxi=max(max(map));
% elseif isochoice == 2
% mini=str2num(get(handles.isomin,'String'));
% maxi=str2num(get(handles.isomax,'String'));
% end
%
% imshow(map, [0 maxi], 'InitialMagnification', 800),
% hold on
% jetcolormap = (colormap('jet'));
% jetcolormap(1,:) = [1, 1, 1];
% colormap(jetcolormap);
% caxis([mini maxi]);
% freezeColors
% hold on
% plot(actcol,actrow,'r+','MarkerSize',20,'LineWidth',3);
% hold off
% title('Please choose centre and press enter');
% [rcol rrow] = getpts(gcf)
% rrow=round(rrow)
% rcol=round(rcol)
% if rrow(1) == 0
%     rrow(1) = 1
% end
%
% if rcol(1) == 0
%     rcol(1) = 1
% end
% hold on
% plot(rcol,rrow,'k+','MarkerSize',20,'LineWidth',3);
% th = 0:pi/50:2*pi;
% xunit = veldistancepix * cos(th) + rcol;
% yunit = veldistancepix * sin(th) + rrow;
% h = plot(xunit, yunit,'k','LineWidth',3);
%
%
% count=0;
% for i =1:numel(xunit)
%
%     time_A=(map(rrow,rcol));
%     time_B=(map(round(yunit(i)),round(xunit(i))));
%     p1=[rcol,rrow];
%     p2=[xunit(i),yunit(i)];
%     tim=time_B-time_A;
%     if tim > 0
%         count=count+1;
%         mat(count,1)=xunit(i);
%         mat(count,2)=yunit(i);
%         mat(count,3)=tim;
%         mat(count,4)=(veldistancepix*0.0001*(str2num(get(handles.pixelsize,'String'))))/tim*1000;
%     end
% end
% [tim_long, ind_long] = max(mat(:,3));
% [tim_short, ind_short] = min(mat(:,3));
%
% p2slow=[xunit(ind_long),yunit(ind_long)];
% p2quick=[xunit(ind_short),yunit(ind_short)];
% hold on
% dp=p2quick-p1;
% q=quiver(p1(1),p1(2),dp(1),dp(2),0);
% q.Color='b';
% q.LineWidth=3;
%
%
% dp=p2slow-p1;
% q=quiver(p1(1),p1(2),dp(1),dp(2),0);
% q.Color='r';
% q.LineWidth=3;
% hold off
%
% tim_avg=mean(mat(:,3));
%
%
% speed_quick=mat(ind_short,4);
% speed_slow=mat(ind_long,4);
% speed_avg=mean(mat(:,4));
% h=msgbox({'Longest Activation Time =';num2str(tim_long);'Shortest Activation Time =';num2str(tim_short);'Average Activation Time =';num2str(tim_avg)})
% h=msgbox({'Slowest Speed =';num2str(speed_slow);'Fastest Speed';num2str(speed_quick);'Average Speed =';num2str(speed_avg)})
%
%     case 'Activation Curve'
% [map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,~,vout,quivers_Xout,quivers_Yout,quivers_vxout,quivers_vyout, onedevcv]...
% = cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'));
%
% tim=act_t;
% tim=tim-min(tim);
% allpts=numel(tim);
% xbins=0:0.01:max(tim);
% tissueact=100*cumsum(hist(tim,xbins))/allpts;
% figure,
% plot(xbins,tissueact,'k')
% xlabel('Time (ms)');
% ylabel('Tissue Activated (%)');
% title('Activation Curve');
% ax=gca;
% axis tight
% ylim([0 100]);
% xx=get(ax,'Xlim');
% yy=get(ax,'Ylim');
% hold on
%
%
% I90 = find(tissueact > 90);
% I50 = find(tissueact > 50);
% I20 = find(tissueact > 20);
%
% I90=I90(1);
% I50=I50(1);
% I20=I20(1);
%
% tim100=max(xbins);
% tim90=I90*0.01;
% tim50=I50*0.01;
% tim20=I20*0.01;
%
% closest90=tissueact(I90);
% closest50=tissueact(I50);
% closest20=tissueact(I20);
%
% ha20 = plot([tim20 tim20], [yy(1) closest20] ,'LineStyle', '--', 'Color', 'g');
% hb20 = plot([xx(1) tim20], [closest20 closest20],'LineStyle', '--', 'Color', 'g');
%
% ha50 = plot([tim50 tim50], [yy(1) closest50] ,'LineStyle', '--', 'Color', 'b');
% hb50 = plot([xx(1) tim50], [closest50 closest50],'LineStyle', '--', 'Color', 'b');
%
% ha90 = plot([tim90 tim90], [yy(1) closest90] ,'LineStyle', '--', 'Color', 'r');
% hb90 = plot([xx(1) tim90], [closest90 closest90],'LineStyle', '--', 'Color', 'r');
%
%
% text(0.6,0.4,['20% activation at ', num2str(tim20),' ms'],'Units','normalized');
% text(0.6,0.3,['50% activation at ', num2str(tim50),' ms'],'Units','normalized');
% text(0.6,0.2,['90% activation at ', num2str(tim90),' ms'],'Units','normalized');
% text(0.6,0.1,['100% activation at ', num2str(tim100),' ms'],'Units','normalized');
% end



function beforeGUI_Callback(hObject, eventdata, handles)
% hObject    handle to beforeGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beforeGUI as text
%        str2double(get(hObject,'String')) returns contents of beforeGUI as a double


% --- Executes during object creation, after setting all properties.
function beforeGUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beforeGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function afterGUI_Callback(hObject, eventdata, handles)
% hObject    handle to afterGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of afterGUI as text
%        str2double(get(hObject,'String')) returns contents of afterGUI as a double


% --- Executes during object creation, after setting all properties.
function afterGUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to afterGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in custfilt.
function custfilt_Callback(hObject, eventdata, handles)
% hObject    handle to custfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filterBuilder


% --- Executes on selection change in velalgo.
function velalgo_Callback(hObject, eventdata, handles)
% hObject    handle to velalgo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns velalgo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velalgo


% --- Executes during object creation, after setting all properties.
function velalgo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velalgo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in squareROI.
function squareROI_Callback(hObject, eventdata, handles)
% hObject    handle to squareROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of squareROI



function squaresize_Callback(hObject, eventdata, handles)
% hObject    handle to squaresize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of squaresize as text
%        str2double(get(hObject,'String')) returns contents of squaresize as a double


% --- Executes during object creation, after setting all properties.
function squaresize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to squaresize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in actopt.
function actopt_Callback(hObject, eventdata, handles)
% hObject    handle to actopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns actopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from actopt


% --- Executes during object creation, after setting all properties.
function actopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function actmin_Callback(hObject, eventdata, handles)
% hObject    handle to actmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,~,vout,quivers_Xout,quivers_Yout,quivers_vxout,quivers_vyout, onedevcv]...
    = cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
    str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));
tim=act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2num(get(handles.actmax,'String'));
actmin=str2num(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);

if actmax < 100
    Imax=Imax(1);
else Imax=max(tim);
end
Imin=Imin(1);


if actmax < 100
    timmax=Imax*0.01;
else timmax=Imax
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
    timdiff=timmax;
end
if get(handles.checkbox8,'Value')==1
    normfac=225/allpts;
    timdiff=timdiff*normfac;
end

set(handles.actquote, 'String', [num2str(timdiff),' ms']);

guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of actmin as text
%        str2double(get(hObject,'String')) returns contents of actmin as a double


% --- Executes during object creation, after setting all properties.
function actmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function actmax_Callback(hObject, eventdata, handles)
% hObject    handle to actmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of actmax as text
%        str2double(get(hObject,'String')) returns contents of actmax as a double

handles = guidata(hObject);
[map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,~,vout,quivers_Xout,quivers_Yout,quivers_vxout,quivers_vyout, onedevcv]...
    = cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
    str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));

tim=act_t;
tim=tim-min(tim);
allpts=numel(tim);
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2num(get(handles.actmax,'String'));
actmin=str2num(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);

if actmax < 100
    Imax=Imax(1);
else Imax=max(tim);
end

Imin=Imin(1);


if actmax < 100
    timmax=Imax*0.01;
else timmax=Imax
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
    timdiff=timmax;
end

if get(handles.checkbox8,'Value')==1
    normfac=225/allpts;
    timdiff=timdiff*normfac;
end
set(handles.actquote, 'String', [num2str(timdiff),' ms']);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function actmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[map,act_x,act_y,act_t,quivers_X,quivers_Y,quivers_vx,quivers_vy,CV,~,vout,quivers_Xout,quivers_Yout,quivers_vxout,quivers_vyout, onedevcv]...
    = cvmap(str2num(get(handles.pixelsize,'String')),str2num(get(handles.framerate,'String')),handles.cvimages,handles.mask,get(handles.velout,'Value'),str2num(get(handles.minvel,'String')),str2num(get(handles.maxvel,'String')),get(handles.velalgo,'Value'),...
    str2num(get(handles.MINt,'String')),str2num(get(handles.MAXt,'String')),str2num(get(handles.winsize,'String')),str2num(get(handles.beforeGUI,'String')),str2num(get(handles.wint,'String')),0,str2num(get(handles.t,'String')),get(handles.tfilt,'Value'),get(handles.usespline,'Value'),str2num(get(handles.splineN,'String')));

tim=act_t;
tim=tim-min(tim);
allpts=numel(tim)
xbins=0:0.01:max(tim);
tissueact=100*cumsum(hist(tim,xbins))/allpts;

actmax=str2num(get(handles.actmax,'String'));
actmin=str2num(get(handles.actmin,'String'));

Imax = find(tissueact > actmax);
Imin = find(tissueact > actmin);

if actmax < 100
    Imax=Imax(1);
else Imax=max(tim);
end
Imin=Imin(1);


if actmax < 100
    timmax=Imax*0.01;
else timmax=Imax
end
timmin=Imin*0.01;
timdiff=timmax-timmin;
if actmin == 0
    timdiff=timmax;
end
if get(handles.checkbox8,'Value')==1
    normfac=225/allpts;
    timdiff=timdiff*normfac;
end

set(handles.actquote, 'String', [num2str(timdiff),' ms']);
% Hint: get(hObject,'Value') returns toggle state of checkbox8



function t_Callback(hObject, eventdata, handles)
% hObject    handle to t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of t as text
%        str2double(get(hObject,'String')) returns contents of t as a double


% --- Executes during object creation, after setting all properties.
function t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apdscale.
function apdscale_Callback(hObject, eventdata, handles)
% hObject    handle to apdscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns apdscale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apdscale


% --- Executes during object creation, after setting all properties.
function apdscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MINt_Callback(hObject, eventdata, handles)
% hObject    handle to MINt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MINt as text
%        str2double(get(hObject,'String')) returns contents of MINt as a double


% --- Executes during object creation, after setting all properties.
function MINt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MINt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MAXt_Callback(hObject, eventdata, handles)
% hObject    handle to MAXt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MAXt as text
%        str2double(get(hObject,'String')) returns contents of MAXt as a double


% --- Executes during object creation, after setting all properties.
function MAXt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MAXt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in actfittimes.
function actfittimes_Callback(hObject, eventdata, handles)
% hObject    handle to actfittimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns actfittimes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from actfittimes


% --- Executes during object creation, after setting all properties.
function actfittimes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actfittimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function winsize_Callback(hObject, eventdata, handles)
% hObject    handle to winsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winsize as text
%        str2double(get(hObject,'String')) returns contents of winsize as a double


% --- Executes during object creation, after setting all properties.
function winsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function taustart_Callback(hObject, eventdata, handles)
% hObject    handle to taustart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
contents = cellstr(get(handles.Mapchoice,'String'));
choice=get(handles.Mapchoice,'Value');
if choice == 8
    Mapchoice_Callback(hObject, eventdata, handles)
end
% Hints: get(hObject,'String') returns contents of taustart as text
%        str2double(get(hObject,'String')) returns contents of taustart as a double


% --- Executes during object creation, after setting all properties.
function taustart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to taustart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r2cut_Callback(hObject, eventdata, handles)
% hObject    handle to r2cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
contents = cellstr(get(handles.Mapchoice,'String'));
choice=get(handles.Mapchoice,'Value');
if choice == 8
    Mapchoice_Callback(hObject, eventdata, handles)
end
% Hints: get(hObject,'String') returns contents of r2cut as text
%        str2double(get(hObject,'String')) returns contents of r2cut as a double


% --- Executes during object creation, after setting all properties.
function r2cut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r2cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function apdblnum_Callback(hObject, eventdata, handles)
% hObject    handle to apdblnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of apdblnum as text
%        str2double(get(hObject,'String')) returns contents of apdblnum as a double


% --- Executes during object creation, after setting all properties.
function apdblnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdblnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apdbl.
function apdbl_Callback(hObject, eventdata, handles)
% hObject    handle to apdbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns apdbl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apdbl


% --- Executes during object creation, after setting all properties.
function apdbl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function taufinish_Callback(hObject, eventdata, handles)
% hObject    handle to taufinish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of taufinish as text
%        str2double(get(hObject,'String')) returns contents of taufinish as a double


% --- Executes during object creation, after setting all properties.
function taufinish_CreateFcn(hObject, eventdata, handles)
% hObject    handle to taufinish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function peakhigh_Callback(hObject, eventdata, handles)
% hObject    handle to peakhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of peakhigh as text
%        str2double(get(hObject,'String')) returns contents of peakhigh as a double


% --- Executes during object creation, after setting all properties.
function peakhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peakhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scal_Callback(hObject, eventdata, handles)
% hObject    handle to scal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scal as text
%        str2double(get(hObject,'String')) returns contents of scal as a double


% --- Executes during object creation, after setting all properties.
function scal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wint_Callback(hObject, eventdata, handles)
% hObject    handle to wint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wint as text
%        str2double(get(hObject,'String')) returns contents of wint as a double


% --- Executes during object creation, after setting all properties.
function wint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in altanal.
function altanal_Callback(hObject, eventdata, handles)
% hObject    handle to altanal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
alternangui


% --- Executes on button press in removef.
function removef_Callback(hObject, eventdata, handles)
% hObject    handle to removef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of removef


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dual2


% --- Executes on selection change in aptime1.
function aptime1_Callback(hObject, eventdata, handles)
% hObject    handle to aptime1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns aptime1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from aptime1


% --- Executes during object creation, after setting all properties.
function aptime1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aptime1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in configure.
function configure_Callback(hObject, eventdata, handles)
% hObject    handle to configure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Construct a questdlg
choice = questdlg('Would you like to load a configuration file or save current settings?', ...
    'Config File', ...
    'New File', 'Load Settings','Load Settings');
% Handle response
switch choice
    case 'New File'
        [fname,PathName] = uiputfile('*.txt')
        filename=[PathName,fname];
        
        % create file for writing too
        fileID = fopen(filename,'w');
        dt=datetime('now'); dt=datestr(dt);
        fprintf(fileID,'Configuration File for use in ElectroMap\r\n');
        fprintf(fileID,['Created: ',dt,'\r\n'])
        fprintf(fileID,['Notes:\r\n\r\n\r\n'])
        fprintf(fileID,['--- DO NOT EDIT VARIABLE NAMES OR REMOVE ! BELOW THIS POINT, {} = units, () = settings in ElectroMap or Notes ---\r\n\r\n'])
        
        %cell arrays for popdown meny strings
        threshopt_string=get(handles.threshopt,'String');
        thershopt_string=threshopt_string{get(handles.threshopt,'Value')};
        sfilt_string=get(handles.sfilt,'String');
        sfilt_string=sfilt_string{get(handles.sfilt,'Value')};
        segchoice_string=get(handles.segchoice,'String');
        segchoice_string=segchoice_string{get(handles.segchoice,'Value')}
        BLopt_string=get(handles.BLopt,'String');
        BLopt_string=BLopt_string{get(handles.BLopt,'Value')}
        tfilt_string=get(handles.tfilt,'String');
        tfilt_string=tfilt_string{get(handles.tfilt,'Value')};
        apdbl_string=get(handles.apdbl,'String');
        apdbl_string=apdbl_string{get(handles.apdbl,'Value')};
        aptime1_string=get(handles.aptime1,'String');
        aptime1_string=aptime1_string{get(handles.aptime1,'Value')};
        velalgo_string=get(handles.velalgo,'String');
        velalgo_string=velalgo_string{get(handles.velalgo,'Value')};
        actfittimes_string=get(handles.actfittimes,'String');
        actfittimes_string=actfittimes_string{get(handles.actfittimes,'Value')};
        velout_string=get(handles.velout,'String');
        velout_string=velout_string{get(handles.velout,'Value')};
        apdout_string=get(handles.apdout,'String');
        apdout_string=apdout_string{get(handles.apdout,'Value')};
        apdscale_string=get(handles.apdscale,'String');
        apdscale_string=apdscale_string{get(handles.apdscale,'Value')};
        % get settings and put into file
        fprintf(fileID,['framerate=',get(handles.framerate,'String'),'! {kHz} !\r\n'])
        fprintf(fileID,['pixelsize=',get(handles.pixelsize,'String'),'! {ms} !!\r\n'])
        fprintf(fileID,['threshopt=',num2str(get(handles.threshopt,'Value')),'! (',thershopt_string,')',' !\r\n'])
        fprintf(fileID,['manthresh=',get(handles.manthresh,'String'),'! {Percent} (Change from automatically generated threshold) !\r\n'])
        fprintf(fileID,['sfilt=',num2str(get(handles.sfilt,'Value')),'! (',sfilt_string,') !\r\n'])
        fprintf(fileID,['sfiltsize=',get(handles.sfiltsize,'String'),'! {Pixels} !\r\n'])
        fprintf(fileID,['segchoice=',num2str(get(handles.segchoice,'Value')),'! (',segchoice_string,') !\r\n'])
        fprintf(fileID,['segsize=',get(handles.segsize,'String'),'! !\r\n'])
        fprintf(fileID,['invertopt=',num2str(get(handles.invertopt,'Value')),'! (1 means invert signal) !\r\n'])
        fprintf(fileID,['BLopt=',num2str(get(handles.BLopt,'Value')),'! (',BLopt_string,') !\r\n'])
        fprintf(fileID,['thlen=',get(handles.thlen,'String'),'! {ms}, (Length of top-hat filter) !\r\n'])
        fprintf(fileID,['tfilt=',num2str(get(handles.tfilt,'Value')),'! (',tfilt_string,') !\r\n'])
        fprintf(fileID,['minpeak=',get(handles.minpeak,'String'),'! {ms} !\r\n'])
        fprintf(fileID,['peakhigh=',get(handles.peakhigh,'String'),'! !\r\n'])
        fprintf(fileID,['minnum=',get(handles.minnum,'String'),'! !\r\n'])
        fprintf(fileID,['minbound=',get(handles.minbound,'String'),'! {ms} !\r\n'])
        fprintf(fileID,['beforeGUI=',get(handles.beforeGUI,'String'),'! {ms} !\r\n'])
        fprintf(fileID,['afterGUI=',get(handles.afterGUI,'String'),'! {ms} !\r\n'])
        fprintf(fileID,['apdbl=',num2str(get(handles.apdbl,'Value')),'! (',apdbl_string,') !\r\n'])
        fprintf(fileID,['apdblnum=',get(handles.apdblnum,'String'),'! {ms} !\r\n'])
        fprintf(fileID,['aptime1=',num2str(get(handles.aptime1,'Value')),'! (',aptime1_string,') !\r\n'])
        fprintf(fileID,['taustart=',get(handles.taustart,'String'),'! {Percent} !\r\n'])
        fprintf(fileID,['taufinish=',get(handles.taufinish,'String'),'! {Percent} !\r\n'])
        fprintf(fileID,['r2cut=',get(handles.r2cut,'String'),'! !\r\n'])
        fprintf(fileID,['velalgo=',num2str(get(handles.velalgo,'Value')),'! (',velalgo_string,') (Activation Measure) !\r\n'])
        fprintf(fileID,['isoopt=',num2str(get(handles.isoopt,'Value')),'! !\r\n'])
        fprintf(fileID,['isomin=',get(handles.isomin,'String'),'! {ms} !\r\n'])
        fprintf(fileID,['isomax=',get(handles.isomax,'String'),'! {ms} !\r\n'])
        fprintf(fileID,['actfittimes=',num2str(get(handles.actfittimes,'Value')),'! (',actfittimes_string,') !\r\n'])
        fprintf(fileID,['MINt=',get(handles.MINt,'String'),'! {ms} (Minimum activation time for multi-vector fit) !\r\n'])
        fprintf(fileID,['MAXt=',get(handles.MAXt,'String'),'! {ms} (Maximum activation time for multi-vector fit) !\r\n'])
        fprintf(fileID,['velout=',num2str(get(handles.velout,'Value')),'! (',velout_string,') (Local Velocity outlier removal) !\r\n'])
        fprintf(fileID,['minvel=',get(handles.minvel,'String'),'! {ms} (Minimum calcualted velocity that is not discarded) !\r\n'])
        fprintf(fileID,['maxvel=',get(handles.maxvel,'String'),'! {ms} (Maximum calculated velocity that is not discarded) !\r\n'])
        fprintf(fileID,['winsize=',get(handles.winsize,'String'),'!{Pixels} (Local window size) !\r\n'])
        fprintf(fileID,['scal=',get(handles.scal,'String'),'! (Size of overlaid velocity vectors) !\r\n'])
        fprintf(fileID,['wint=',get(handles.wint,'String'),'! (Maximum time diffrence allowed in local window fit) !\r\n'])
        fprintf(fileID,['apdout=',num2str(get(handles.apdout,'Value')),'! (',apdout_string,') (APD/CaD outlier removal) !\r\n'])
        fprintf(fileID,['apdscale=',num2str(get(handles.apdscale,'Value')),'! (',apdscale_string,') !\r\n'])
        fprintf(fileID,['cmin=',get(handles.cmin,'String'),'! {ms} (manual colour map minimum)!\r\n'])
        fprintf(fileID,['cmax=',get(handles.cmax,'String'),'! {ms} (manual colour map maximum)!\r\n'])
        fprintf(fileID,['t=',get(handles.t,'String'),'! {Percent} (APD/CaD)!\r\n'])
        fprintf(fileID,['checkbox8=',num2str(get(handles.checkbox8,'Value')),'! (1 means normalised to {ms/mm2}, 0 absoulte in {ms})!\r\n'])
        fprintf(fileID,['actmin=',get(handles.actmin,'String'),'! {Percent}!\r\n'])
        fprintf(fileID,['actmax=',get(handles.actmax,'String'),'! {Percent}!\r\n'])
        fprintf(fileID,['binnumber=',get(handles.binnumber,'String'),'!!\r\n'])
        
        % close file
        fclose(fileID)
        
    case 'Load Settings'
        [fname,PathName] = uigetfile('*.txt')
        filename=[PathName,fname]
        %% open file for reading
        fid=fopen(filename,'r','b')
        fstr=fread(fid,'int8=>char')';
        fclose(fid);
        
        %% Update GUI
        set(handles.framerate,'String',varEM(fstr,'framerate',0))
        set(handles.pixelsize,'String',varEM(fstr,'pixelsize',0))
        set(handles.threshopt,'Value',varEM(fstr,'threshopt',1))
        set(handles.manthresh,'String',varEM(fstr,'manthresh',0))
        set(handles.sfilt,'Value',varEM(fstr,'sfilt',1))
        set(handles.sfiltsize,'String',varEM(fstr,'sfiltsize',0))
        set(handles.segchoice,'Value',varEM(fstr,'segchoice',1))
        set(handles.segsize,'String',varEM(fstr,'segsize',0))
        set(handles.invertopt,'Value',varEM(fstr,'invertopt',1))
        set(handles.BLopt,'Value',varEM(fstr,'BLopt',1))
        set(handles.tfilt,'Value',varEM(fstr,'tfilt',1))
        set(handles.minpeak,'String',varEM(fstr,'minpeak',0))
        set(handles.peakhigh,'String',varEM(fstr,'peakhigh',0))
        set(handles.minnum,'String',varEM(fstr,'minnum',0))
        set(handles.minbound,'String',varEM(fstr,'minbound',0))
        set(handles.beforeGUI,'String',varEM(fstr,'beforeGUI',0))
        set(handles.afterGUI,'String',varEM(fstr,'afterGUI',0))
        set(handles.apdbl,'Value',varEM(fstr,'apdbl',1))
        set(handles.apdblnum,'String',varEM(fstr,'apdblnum',0))
        set(handles.aptime1,'Value',varEM(fstr,'aptime1',1))
        set(handles.taustart,'String',varEM(fstr,'taustart',0))
        set(handles.taufinish,'String',varEM(fstr,'taufinish',0))
        set(handles.r2cut,'String',varEM(fstr,'r2cut',0))
        set(handles.velalgo,'Value',varEM(fstr,'velalgo',1))
        set(handles.isoopt,'Value',varEM(fstr,'isoopt',1))
        set(handles.isomin,'String',varEM(fstr,'isomin',0))
        set(handles.isomax,'String',varEM(fstr,'isomax',0))
        set(handles.actfittimes,'Value',varEM(fstr,'actfittimes',1))
        set(handles.MINt,'String',varEM(fstr,'MINt',0))
        set(handles.MAXt,'String',varEM(fstr,'MAXt',0))
        set(handles.velout,'Value',varEM(fstr,'velout',1))
        set(handles.minvel,'String',varEM(fstr,'minvel',0))
        set(handles.maxvel,'String',varEM(fstr,'maxvel',0))
        set(handles.winsize,'String',varEM(fstr,'winsize',0))
        set(handles.scal,'String',varEM(fstr,'scal',0))
        set(handles.wint,'String',varEM(fstr,'wint',0))
        set(handles.apdout,'Value',varEM(fstr,'apdout',1))
        set(handles.apdscale,'Value',varEM(fstr,'apdscale',1))
        set(handles.cmin,'String',varEM(fstr,'cmin',0))
        set(handles.cmax,'String',varEM(fstr,'cmax',0))
        set(handles.t,'String',varEM(fstr,'t',0))
        set(handles.checkbox8,'Value',varEM(fstr,'checkbox8',1))
        set(handles.actmin,'String',varEM(fstr,'actmin',0))
        set(handles.actmax,'String',varEM(fstr,'actmax',0))
        set(handles.binnumber,'String',varEM(fstr,'binnumber',0))
        set(handles.thlen,'String',varEM(fstr,'thlen',0))
end



function thlen_Callback(hObject, eventdata, handles)
% hObject    handle to thlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thlen as text
%        str2double(get(hObject,'String')) returns contents of thlen as a double


% --- Executes during object creation, after setting all properties.
function thlen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in colmap.
function colmap_Callback(hObject, eventdata, handles)
% hObject    handle to colmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns colmap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colmap


% --- Executes during object creation, after setting all properties.
function colmap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in drawcon.
function drawcon_Callback(hObject, eventdata, handles)
% hObject    handle to drawcon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
Mapchoice_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of drawcon




% --- Executes on button press in roibutton.
function roibutton_Callback(hObject, eventdata, handles)
% hObject    handle to roibutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Construct a questdlg with three options
handles = guidata(hObject);
handles
choice = questdlg('Save Current ROI or load previous?', ...
    'ROI', ...
    'Save ROI','Load ROI','Save ROI');
% Handle response
switch choice
    case 'Save ROI'
        [filename,pathname] = uiputfile({'*.txt'}, 'Save ROI in text file');
        [~,~,ext] = fileparts(filename);
        file=[pathname,filename];
        savemask=handles.mask;
        figure,
        imshow(savemask,[])
        dlmwrite(file,savemask);
        lmask=maskload(file)
        figure,
        imshow(lmask,[])
        size(savemask)
        size(lmask)
    case 'Load ROI'
        [filename,pathname] = uigetfile('*.txt','Select the ROI File');
        file=[pathname,filename];
        handles.loadedmask=maskload(file)
        handles.herefromroiload=1;
        guidata(hObject, handles);
        pushload_Callback(hObject, eventdata, handles)
        handles = guidata(hObject); %%update handles after doing OMimload
        handles.herefromroiload=0;
        guidata(hObject, handles);
end
guidata(hObject, handles);


% --- Executes on button press in rawvid.
function rawvid_Callback(hObject, eventdata, handles)
% hObject    handle to rawvid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
axes(handles.imageaxes)
background = repmat(handles.im, [1, 1, 3]);
[rows, cols]=size(handles.im);
images=handles.images;
mask=handles.mask;
images=imcomplement(images);
im=handles.im;
im=double(im);
im=im-min(min(im));
im=im./max(max(im));
im=im*65535;
im=uint16(im);

savegif=1;
if savegif == 1
    [a,b]=uiputfile('*.gif');
    filename=[b,a];
end
prompt = {'Fluorescence threshold (0-1):','Video Start (s):','Video End (s)','Normalise? (0=no, 1=yes)'};
dlg_title = 'Raw Video Options';
num_lines = 1;
exposure=1/str2num(get(handles.framerate,'String'));
defaultans = {'0.2','0',num2str(size(images,3)/1000*exposure),'1'};
opts = inputdlg(prompt,dlg_title,num_lines,defaultans);
flthresh=str2num(opts{1});
istart=str2num(opts{2});

iend=str2num(opts{3});
normF=str2num(opts{4});
%change is to frame #
istart=round((istart/exposure)*1000);
iend=round((iend/exposure)*1000);
if istart == 0
    istart = 1
end

if iend > size(images,3)
    iend = size(images,3)
end

for r=1:rows
    for c=1:cols
        sig=[];
        sig=squeeze(images(r,c,:));
        sig=sig-min(sig);
        sig=double(sig);
        if normF == 1
            sig=(sig./max(sig))*65535;
        end
        sig=uint16(sig);
        images(r,c,:)=sig;
    end
end
wb=waitbar(0,'Saving Raw Video');
images=double(images);
mask=double(mask);
maxval=max(max(max(images)));
background = repmat(im, [1, 1, 3]);
for i =istart:iend
    waitbar(i/(iend-istart),wb,'Saving Raw Video');
    combinedImage = background;
    foreground=images(:,:,i).*mask;
    foreground=foreground./maxval;
    foreground(foreground < flthresh) = 0;
    %pause(5)
    foregroundColourised = colouriseData(foreground, 'j',flthresh,1);
    c1 = combinedImage(:, :, 1);
    c2 = combinedImage(:, :, 2);
    c3 = combinedImage(:, :, 3);
    
    f1 = foregroundColourised(:, :, 1);
    f2 = foregroundColourised(:, :, 2);
    f3 = foregroundColourised(:, :, 3);
    
    c1(sum(foreground, 3) ~= 0) = f1(sum(foreground, 3) ~= 0);
    c2(sum(foreground, 3) ~= 0) = f2(sum(foreground, 3) ~= 0);
    c3(sum(foreground, 3) ~= 0) = f3(sum(foreground, 3) ~= 0);
    
    combinedImage(:, :, 1) = c1;
    combinedImage(:, :, 2) = c2;
    combinedImage(:, :, 3) = c3;
    
    hold on
    
    axis image;
    axis off;
    
    if savegif == 1
        delay=0.01;
        [imind,cm] = rgb2ind(combinedImage,256);
        % Write to the GIF File
        if i == istart
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay);
        end
    end
end
delete(wb)

% --- Executes on button press in resegment.
function resegment_Callback(hObject, eventdata, handles)
% hObject    handle to resegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.herefromsegmentpush=1;
guidata(hObject, handles);
pushbutton12_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
handles.herefromsegmentpush = 0;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
guidata(hObject, handles);

% --- Executes on button press in B2B.
function B2B_Callback(hObject, eventdata, handles)
% hObject    handle to B2B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.herefromsegmentpush=1;
set(handles.segsize,'String',1);
set(handles.segchoice,'Value',2);
guidata(hObject, handles);
pushbutton12_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
handles.herefromsegmentpush = 0;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
guidata(hObject, handles);


% --- Executes on selection change in winopt.
function winopt_Callback(hObject, eventdata, handles)
% hObject    handle to winopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns winopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from winopt


% --- Executes during object creation, after setting all properties.
function winopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segsignal.
function segsignal_Callback(hObject, eventdata, handles)
% hObject    handle to segsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.herefromsegmentpush=1;
guidata(hObject, handles);
pushbutton12_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
% Hints: contents = cellstr(get(hObject,'String')) returns segsignal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from segsignal


% --- Executes during object creation, after setting all properties.
function segsignal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sfiltsigma_Callback(hObject, eventdata, handles)
% hObject    handle to sfiltsigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sfiltsigma as text
%        str2double(get(hObject,'String')) returns contents of sfiltsigma as a double


% --- Executes during object creation, after setting all properties.
function sfiltsigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfiltsigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in usespline.
function usespline_Callback(hObject, eventdata, handles)
% hObject    handle to usespline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns usespline contents as cell array
%        contents{get(hObject,'Value')} returns selected item from usespline


% --- Executes during object creation, after setting all properties.
function usespline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usespline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function splineN_Callback(hObject, eventdata, handles)
% hObject    handle to splineN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of splineN as text
%        str2double(get(hObject,'String')) returns contents of splineN as a double


% --- Executes during object creation, after setting all properties.
function splineN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to splineN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function freqmapopt_Callback(hObject, eventdata, handles)
% hObject    handle to freqmapopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles = guidata(hObject); 
    prompt = {'Minimum Frequency (Hz):','Maximum Frequency (Hz):','Frequency Bin Size (Hz)','Window? 0 = no, 1 = hann'};
    dims = [1 35];
    definput = {num2str(handles.fmin),num2str(handles.fmax),num2str(handles.fbin),num2str(handles.dfwin)};
    answer = inputdlg(prompt,'Frequnecy Mapping Options',dims,definput)
    handles.fmin=str2num(answer{1});
    handles.fmax=str2num(answer{2});
    handles.fbin=str2num(answer{3});
    handles.dfwin=str2num(answer{4});
    guidata(hObject, handles);
    if get(handles.Mapchoice,'Value') == 5;
        Mapchoice_Callback(hObject, eventdata, handles)
    end
    guidata(hObject, handles);

    


% --------------------------------------------------------------------
function ROInum_Callback(hObject, eventdata, handles)
% hObject    handle to ROInum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
    prompt = {'Number of ROIs:','Remove overlapping pixels? (0=Yes) (1=No)'};
    dims = [1 35];
    definput = {num2str(handles.roinum),num2str(handles.roisum)};
    answer = inputdlg(prompt,'ROI options',dims,definput)
    handles.roinum=str2num(answer{1});
    handles.roisum=str2num(answer{2});
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function coljet_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
set(handles.colmap,'Value',1);
 Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);



% --------------------------------------------------------------------
function colhsv_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
set(handles.colmap,'Value',2);
 Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --------------------------------------------------------------------
function colhot_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
set(handles.colmap,'Value',3);
 Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colcool_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
set(handles.colmap,'Value',4);
 Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colparula_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
set(handles.colmap,'Value',5);
 Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --------------------------------------------------------------------
function colspring_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
set(handles.colmap,'Value',6);
 Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colsummer_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
set(handles.colmap,'Value',7);
 Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colautumn_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
set(handles.colmap,'Value',8);
 Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function colwinter_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
set(handles.colmap,'Value',9);
 Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function bgblack_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
handles.bgcol='k';
handles.bgon=1;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --------------------------------------------------------------------
function bgwhite_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
handles.bgcol='w';
handles.bgon=1;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --------------------------------------------------------------------
function bgtran_Callback(hObject, eventdata, handles)
handles = guidata(hObject); 
handles.bgon=0;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function snrcalc_Callback(hObject, eventdata, handles)
% hObject    handle to snrcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ttpset_Callback(hObject, eventdata, handles)
% hObject    handle to ttpset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles = guidata(hObject); 
    prompt = {'Start Point (%):','End Point (%):'};
    dims = [1 35];
    definput = {num2str(handles.ttpstart),num2str(handles.ttpend)};
    answer = inputdlg(prompt,'Frequnecy Mapping Options',dims,definput);
    handles.ttpstart=str2num(answer{1});
    handles.ttpend=str2num(answer{2});
    guidata(hObject, handles);
    if get(handles.Mapchoice,'Value') == 7;
        Mapchoice_Callback(hObject, eventdata, handles)
    end
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function connnnnnnnnnn_Callback(hObject, eventdata, handles)
% hObject    handle to connnnnnnnnnn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function conoff_Callback(hObject, eventdata, handles)
% hObject    handle to conoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.drawcon=0;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function conon_Callback(hObject, eventdata, handles)
% hObject    handle to conon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.drawcon=1;
    prompt = {'Contour spacing (map units):'};
    dims = [1 35];
    if isempty(handles.conbon) == 1
        handles.framerate
        conbount=1/str2num(get(handles.framerate,'String'));
        handles.conbon=num2str(conbount)
    end
    definput = {num2str(handles.conbon)};
    answer = inputdlg(prompt,'Contour Setting',dims,definput);
    handles.conbon=answer{1};
    guidata(hObject, handles);
    Mapchoice_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function nomedifilt_Callback(hObject, eventdata, handles)
% hObject    handle to nomedifilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.medifilt=0;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function yesmedifilt_Callback(hObject, eventdata, handles)
% hObject    handle to yesmedifilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.medifilt=1;
guidata(hObject, handles);
Mapchoice_Callback(hObject, eventdata, handles)
