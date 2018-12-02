function varargout = dual2(varargin)
% DUAL2 MATLAB code for dual2.fig
%      DUAL2, by itself, creates a new DUAL2 or raises the existing
%      singleton*.
%
%      H = DUAL2 returns the handle to a new DUAL2 or the handle to
%      the existing singleton*.
%
%      DUAL2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DUAL2.M with the given input arguments.
%
%      DUAL2('Property','Value',...) creates a new DUAL2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dual2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dual2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dual2

% Last Modified by GUIDE v2.5 09-Aug-2018 18:52:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dual2_OpeningFcn, ...
                   'gui_OutputFcn',  @dual2_OutputFcn, ...
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


% --- Executes just before dual2 is made visible.
function dual2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dual2 (see VARARGIN)
h = findobj('Tag','ElectroMap');
%% get intial data from ElectroMap and set other stuff
g1data = guidata(h);
handles.threshop=get(g1data.threshopt,'Value');
handles.threshman=(get(g1data.manthresh,'Value'));
quinnieopt=0;
handles.rect=[];
handles.row=[];
set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
camopt=0;
handles.framerate=str2num(get(g1data.framerate,'String'));
handles.pixelsize=str2num(get(g1data.pixelsize,'String'));
exposure=1/handles.framerate;
handles.averagestime=g1data.averagestime;
%Baseline
BLopt=(get(g1data.BLopt,'Value'));
handles.q2locs=g1data.q2locs;
handles.dura=str2num(get(g1data.t,'String'));

%% Send 1st file info
handles.fname1=g1data.fname;
handles.frame1=g1data.frame1;
handles.boundaries=g1data.boundaries;
handles.images1=g1data.images;
handles.averages1=g1data.averages;
handles.section=g1data.section;
handles.num_images=g1data.num_images;
handles.mask=g1data.mask;
handles.im1=g1data.im;
handles.I1=g1data.I;
tfilt=get(g1data.tfilt,'Value');
sfilt=get(g1data.sfilt,'Value');
sfiltsize=str2num(get(g1data.sfiltsize,'String'));

%% Setup inputs from ElectroMap

%assume signal 2 opposite to signal 1

inversion=get(g1data.invertopt,'Value')-1;
inversion=abs(inversion);
set(handles.invertopt2,'Value',inversion);

set(handles.cmin1,'String',get(g1data.cmin,'String'));set(handles.cmax1,'String',get(g1data.cmax,'String'));
set(handles.cmin2,'String',get(g1data.cmin,'String'));set(handles.cmax2,'String',get(g1data.cmax,'String'));
set(handles.cmind,'String',get(g1data.cmin,'String'));set(handles.cmaxd,'String',get(g1data.cmax,'String'));

set(handles.conbound1,'String',get(g1data.conbound,'String'));
set(handles.conbound2,'String',get(g1data.conbound,'String'));
set(handles.conboundd,'String',get(g1data.conbound,'String'));

set(handles.colmap,'Value',get(g1data.colmap,'Value'));

set(handles.beforeGUI,'String',get(g1data.beforeGUI,'String'));
set(handles.afterGUI,'String',get(g1data.afterGUI,'String'));

handles.before=round(str2num(get(handles.beforeGUI,'String'))/exposure);
handles.after=round(str2num(get(handles.afterGUI,'String'))/exposure);
%% get 2nd image file
% Construct a questdlg with three options
choice = questdlg('Please choose 2nd File. Does this image require cropping to match size?', ...
	'Load Crop', ...
	'Yes and Process','Yes - with new processing settings','No','Yes');
% Handle response
switch choice
    case 'Yes and Process'
        cropchoice = 1;
        loadprocess=1;

    case 'Yes - with new processing settings'
        cropchoice = 1;
        loadprocess=0;
    case 'No'
        cropchoice = 0;
end
[FileName,PathName] = uigetfile('*.tif;*.TIF;*.tiff;*.TIFF;');
fname=[PathName,FileName];
handles.fname2=fname;
[num_images2,handles.newrect2,mask2,handles.im2,handles.I2,boundaries2,handles.camopt2,handles.frame12,handles.fluoim2] = OMimload(fname,cropchoice,quinnieopt,g1data.threshop,g1data.threshman,handles.rect,get(handles.invertopt2,'Value'),camopt,get(g1data.imagedisp,'Value'));
handles.cropchoice=cropchoice;
%% display images
imchoice=get(handles.imchoice, 'Value');
if imchoice == 1
    axes(handles.map1)
    imshow(handles.frame1,[],'InitialMagnification', 400)
    colormap('gray')
    freezeColors
    hold on
    for i=1:size(handles.boundaries,1)
        plot(handles.boundaries{i}(:,2),handles.boundaries{i}(:,1),'b','LineWidth',2);
    end
    hold off
            axes(handles.cb1);
    cla reset
    hcb=colorbar;
    hcb.Location = 'southoutside'
    ax = gca;
    %axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    axis off
    
    axes(handles.map2)
    imshow(handles.frame12,[],'InitialMagnification', 400)
    colormap('gray')
    freezeColors
    hold on
    for i=1:size(handles.boundaries,1)
        plot(handles.boundaries{i}(:,2),handles.boundaries{i}(:,1),'r','LineWidth',2);
    end
    hold off
    
            axes(handles.cb2);
    cla reset
    hcb=colorbar;
    hcb.Location = 'southoutside'
    ax = gca;
    axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    axis off
end

%% Process 2nd images
if loadprocess == 1
    [handles.preimages2,handles.images2,handles.averages2] = OMimprocess(fname,handles.im2,handles.newrect2,handles.num_images,handles.cropchoice,handles.mask,sfilt,sfiltsize,inversion,tfilt,g1data.frameremove,0,str2num(get(g1data.sfiltsigma,'String')));
    %% Baseline Drift Correction

if BLopt == 1 || BLopt == 4
    th_len=str2num(get(g1data.thlen,'String'));
    th_len=(th_len)/str2num(get(g1data.framerate,'String'));
    th_len=round(th_len);
    se = strel('line', th_len, 0.5);
    BLAV = imopen(handles.averages2, se);
    %figure, plot(BL)
end

if BLopt == 2 || BLopt == 5
    [p,s,mu]=polyfit(1:length(handles.averages2),handles.averages2,4);
    BLAV=polyval(p,1:length(handles.averages2),[],mu);
end

if BLopt == 3 || BLopt == 6
    [p,s,mu]=polyfit(1:length(handles.averages2),handles.averages2,11);
    BLAV=polyval(p,1:length(handles.averages2),[],mu);
end

if BLopt == 7
    BLAV=min(handles.averages2);
end

handles.averages2 = (handles.averages2-BLAV); %Baseline subtraction
%% Remove baseline from each pixel

if BLopt == 4 || BLopt == 5 || BLopt == 6
    for t = 1:size(handles.images2,3)
        %oldsignal(dura)=images(25,50,dura);
        handles.images2(:,:,t)=handles.images2(:,:,t)-BLAV(t);
    end
    for frame = 1:size(handles.images2,3)
        %    newsignal(frame)=images(25,50,frame);
    end
    
end
% if BLopt == 7
% for dura = 1:size(images,3)
%     images(:,:,dura)=images(:,:,dura)-min(min(min(images(:,:,:))));
% end
% end

wb=waitbar(0.5,'Removing Baseline');
revimages2(:,:,:)=imcomplement(handles.images2(:,:,:));

if BLopt == 1 || BLopt == 2 || BLopt == 3
    for row=1:size(handles.images2,1) %%MY BL REMOVAL
        for col=1:size(handles.images2,2)
            for frame = 1:size(handles.images2,3)
                signal(frame)=handles.images2(row,col,frame);
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
            
            
            for frame = 1:size(images,3)
                %oldsignal(frame)=images(25,50,frame);
                handles.images2(row,col,frame)=handles.images2(row,col,frame)+BL(frame);
            end
            handles.images2(row,col,:)=handles.images2(row,col,:) - min(handles.images2(row,col,:)); %make all mins zero
        end
        
    end
    
end
end

%% Signal Axes
axes(handles.avsigs)
hold on
zoom xon
handles.normav1=handles.averages1/max(handles.averages1);
handles.normav2=handles.averages2/max(handles.averages2);
plot(handles.averagestime,handles.normav1,'b')
plot(handles.averagestime,handles.normav2,'r')
xlabel('time(ms)')
set(handles.listbox2,'String',handles.section);
set(handles.listbox2,'Value',get(g1data.listbox2,'Value'));
section_choice=get(handles.listbox2,'Value');
A=(handles.q2locs(section_choice,:));
tstart=min(A(A>0))-handles.before;
if tstart == 0
    tstart = 1;
end
tend=max(A)+handles.after;
if tend > length(handles.averagestime)
    tend = length(handles.averagestime);
end
line([tstart*exposure,tstart*exposure],[0,1],'LineWidth',3,'Color','k')
line([tend*exposure,tend*exposure],[0,1],'LineWidth',3,'Color','k')
% Choose default command line output for dual2
handles.output = hObject;

%% ensemble averaging 1
handles.bframes=ceil(handles.before);
handles.aframes=ceil(handles.after);
APtime = handles.bframes+handles.aframes;
[~,~,num]=size(handles.images1(:,:,:))
size(handles.im1,1)
size(handles.im1,2)
APtime
if APtime <= num
overlay = zeros(size(handles.im1,1), size(handles.im1,2), APtime);
else
overlay = zeros(size(handles.im1,1), size(handles.im1,2), num);
end

m=handles.q2locs(section_choice,:);
f=m(m~=0);
peaks = (length(f));

if f(1) <= handles.bframes
    startloc =2
else startloc =1
end
f;
locRange = startloc:numel(f);
if length(handles.q2locs) >1
for i =1:length(locRange)
    if f(locRange(i)) + handles.aframes < length(handles.images1(1,1,:))
        locRange2(i)=locRange(i);
    end
end
locRange=locRange2;
end

% fill matrix
if length(handles.q2locs) > 1 %only need to overlay if more than 1 peak
    for x = -handles.bframes:handles.aframes
        %if f(locRange)+after < length(handles.images(1,1,:))
        overlay1(:,:,x+handles.bframes+1) = sum(handles.images1(:,:,f(locRange)+x),3)./numel(f);
        overlay1(:,:,x+handles.bframes+1) = overlay1(:,:,x+handles.bframes+1).*double(handles.mask);
        %end
    end
end

if length(handles.q2locs) > 1 %only need to overlay if more than 1 peak
    for x = -handles.bframes:handles.aframes
        %if f(locRange)+after < length(handles.images(1,1,:))
        overlay2(:,:,x+handles.bframes+1) = sum(handles.images2(:,:,f(locRange)+x),3)./numel(f);
        overlay2(:,:,x+handles.bframes+1) = overlay2(:,:,x+handles.bframes+1).*double(handles.mask);
        %end
    end
end


if length(handles.q2locs) == 1 %1 beat
    for x = 1:num
        overlay(:,:,x) = (handles.images(:,:,x));
        overlay(:,:,x) = overlay(:,:,x).*double(handles.mask);
    end
end
handles.cvimages1=overlay1;
handles.cvimages2=overlay2;

    minI1 = min(overlay1(:));
    maxI1 = max(overlay1(:));
    averageBeat1 = overlay1 - minI1;
    averageBeat1 = (2^16-1)*averageBeat1./(maxI1);
    %make 16 bit
    handles.averageBeat1 = uint16(averageBeat1);
    
        minI2 = min(overlay2(:));
    maxI2 = max(overlay2(:));
    averageBeat2 = overlay2 - minI2;
    averageBeat2 = (2^16-1)*averageBeat2./(maxI2);
    %make 16 bit
    handles.averageBeat2 = uint16(averageBeat2);
%% Activation Maps
[handles.actmap1]=activationmap(handles.pixelsize,handles.framerate,handles.cvimages1,handles.mask,get(handles.velalgo,'Value'),str2num(get(handles.beforeGUI,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2num(get(g1data.splineN,'String')));
[handles.actmap2]=activationmap(handles.pixelsize,handles.framerate,handles.cvimages2,handles.mask,get(handles.velalgo,'Value'),str2num(get(handles.beforeGUI,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2num(get(g1data.splineN,'String')));
%% APD/CaD maps

[handles.apmap1,meann,alll,onedev,vari,SE,handles.mapR1]=mapsbabydual(get(handles.velalgo,'Value'),str2num(get(g1data.framerate,'String')),str2num(get(handles.dura2,'String')),handles.I1,handles.images1,handles.averageBeat1,g1data.outlier,str2num(get(handles.cmin1,'String')),str2num(get(handles.cmax1,'String')),get(g1data.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')),get(g1data.medifilt,'Value'));
[handles.apmap2,meann,alll,onedev,vari,SE,handles.mapR2]=mapsbabydual(get(handles.velalgo,'Value'),str2num(get(g1data.framerate,'String')),str2num(get(handles.dura2,'String')),handles.I1,handles.images2,handles.averageBeat2,g1data.outlier,str2num(get(handles.cmin2,'String')),str2num(get(handles.cmax2,'String')),get(g1data.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')),get(g1data.medifilt,'Value'));

%% Dual Maps 
if get(handles.diffopt,'Value') == 1
handles.acttime=handles.actmap1-handles.actmap2;
handles.decaytime=handles.mapR1-handles.mapR2;
handles.duration=handles.apmap1-handles.apmap2;
end

if get(handles.diffopt,'Value') == 2
handles.acttime=handles.actmap2-handles.actmap1;
handles.decaytime=handles.mapR2-handles.mapR1;
handles.duration=handles.apmap2-handles.apmap1;
end
%% Get rid of 0-0 values in dual maps
[rows cols]=size(handles.actmap1)
for r=1:rows
    for c=1:cols
        if handles.actmap1(r,c)==0 || handles.actmap2(r,c)==0
          handles.acttime(r,c)=NaN;
        end
                if handles.mapR1(r,c)==0 || handles.mapR2(r,c)==0
          handles.decaytime(r,c)=NaN;
                end
                if handles.apmap1(r,c)==0 || handles.apmap2(r,c)==0
          handles.duration(r,c)=NaN;
                end 
    end
end
zoom off
% Update handles structure
act_times=handles.acttime(isnan(handles.acttime)==0);
dec_times=handles.decaytime(isnan(handles.decaytime)==0);
dur_times=handles.duration(isnan(handles.duration)==0);

pact=prctile(act_times,[5,50,95]);
pdec=prctile(dec_times,[5,50,95]);
pdur=prctile(dur_times,[5,50,95]);


% Results table
handles.rdata=zeros(3,4);
handles.rdata(1,1)=mean(act_times);handles.rdata(1,2)=std(act_times);handles.rdata(1,3)=std(act_times)/numel(act_times);handles.rdata(1,4)=((pact(3)-pact(1))/pact(2));
handles.rdata(2,1)=mean(dec_times);handles.rdata(2,2)=std(dec_times);handles.rdata(2,3)=std(dec_times)/numel(dec_times);handles.rdata(2,4)=((pdec(3)-pdec(1))/pdec(2));
handles.rdata(3,1)=mean(dur_times);handles.rdata(3,2)=std(dur_times);handles.rdata(3,3)=std(dur_times)/numel(dur_times);handles.rdata(3,4)=((pdur(3)-pdur(1))/pdur(2));
set(handles.rtable,'Data',handles.rdata);


guidata(hObject, handles);


% UIWAIT makes dual2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dual2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
section_choice=get(handles.listbox2,'Value');
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
handles.framerate=str2num(get(g1data.framerate,'String'));
handles.pixelsize=str2num(get(g1data.pixelsize,'String'));
axes(handles.avsigs)
cla
hold on
zoom xon
handles.normav1=handles.averages1/max(handles.averages1);
handles.normav2=handles.averages2/max(handles.averages2);
plot(handles.averagestime,handles.normav1,'b')
plot(handles.averagestime,handles.normav2,'r')
xlabel('time(ms)')

exposure=1/handles.framerate;
A=(handles.q2locs(section_choice,:));
tstart=min(A(A>0))-handles.before;
if tstart == 0
    tstart = 1;
end
tend=max(A)+handles.after;
if tend > length(handles.averagestime)
    tend = length(handles.averagestime);
end
line([tstart*exposure,tstart*exposure],[0,1],'LineWidth',3,'Color','k')
line([tend*exposure,tend*exposure],[0,1],'LineWidth',3,'Color','k')
% Choose default command line output for dual2
handles.output = hObject;

%% ensemble averaging 1
handles.bframes=ceil(handles.before);
handles.aframes=ceil(handles.after);
APtime = handles.bframes+handles.aframes;
[~,~,num]=size(handles.images1(:,:,:))
size(handles.im1,1)
size(handles.im1,2)
APtime
if APtime <= num
overlay = zeros(size(handles.im1,1), size(handles.im1,2), APtime);
else
overlay = zeros(size(handles.im1,1), size(handles.im1,2), num);
end

m=handles.q2locs(section_choice,:);
f=m(m~=0);
peaks = (length(f));

if f(1) <= handles.bframes
    startloc =2
else startloc =1
end
f;
locRange = startloc:numel(f);
if length(handles.q2locs) >1
for i =1:length(locRange)
    if f(locRange(i)) + handles.aframes < length(handles.images1(1,1,:))
        locRange2(i)=locRange(i);
    end
end
locRange=locRange2;
end

% fill matrix
if length(handles.q2locs) > 1 %only need to overlay if more than 1 peak
    for x = -handles.bframes:handles.aframes
        %if f(locRange)+after < length(handles.images(1,1,:))
        overlay1(:,:,x+handles.bframes+1) = sum(handles.images1(:,:,f(locRange)+x),3)./numel(f);
        overlay1(:,:,x+handles.bframes+1) = overlay1(:,:,x+handles.bframes+1).*double(handles.mask);
        %end
    end
end

if length(handles.q2locs) > 1 %only need to overlay if more than 1 peak
    for x = -handles.bframes:handles.aframes
        %if f(locRange)+after < length(handles.images(1,1,:))
        overlay2(:,:,x+handles.bframes+1) = sum(handles.images2(:,:,f(locRange)+x),3)./numel(f);
        overlay2(:,:,x+handles.bframes+1) = overlay2(:,:,x+handles.bframes+1).*double(handles.mask);
        %end
    end
end


if length(handles.q2locs) == 1 %1 beat
    for x = 1:num
        overlay(:,:,x) = (handles.images(:,:,x));
        overlay(:,:,x) = overlay(:,:,x).*double(handles.mask);
    end
end
handles.cvimages1=overlay1;
handles.cvimages2=overlay2;

    minI1 = min(overlay1(:));
    maxI1 = max(overlay1(:));
    averageBeat1 = overlay1 - minI1;
    averageBeat1 = (2^16-1)*averageBeat1./(maxI1);
    %make 16 bit
    handles.averageBeat1 = uint16(averageBeat1);
    
        minI2 = min(overlay2(:));
    maxI2 = max(overlay2(:));
    averageBeat2 = overlay2 - minI2;
    averageBeat2 = (2^16-1)*averageBeat2./(maxI2);
    %make 16 bit
    handles.averageBeat2 = uint16(averageBeat2);
%% Activation Maps
[handles.actmap1]=activationmap(handles.pixelsize,handles.framerate,handles.cvimages1,handles.mask,get(handles.velalgo,'Value'),str2num(get(handles.beforeGUI,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2num(get(g1data.splineN,'String')));
[handles.actmap2]=activationmap(handles.pixelsize,handles.framerate,handles.cvimages2,handles.mask,get(handles.velalgo,'Value'),str2num(get(handles.beforeGUI,'String')),get(g1data.tfilt,'Value'),get(g1data.usespline,'Value'),str2num(get(g1data.splineN,'String')));
%% APD/CaD maps

[handles.apmap1,meann,alll,onedev,vari,SE,handles.mapR1]=mapsbabydual(get(handles.velalgo,'Value'),str2num(get(g1data.framerate,'String')),str2num(get(handles.dura2,'String')),handles.I1,handles.images1,handles.averageBeat1,g1data.outlier,str2num(get(handles.cmin1,'String')),str2num(get(handles.cmax1,'String')),get(g1data.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')),get(g1data.medifilt,'Value'));
[handles.apmap2,meann,alll,onedev,vari,SE,handles.mapR2]=mapsbabydual(get(handles.velalgo,'Value'),str2num(get(g1data.framerate,'String')),str2num(get(handles.dura2,'String')),handles.I1,handles.images2,handles.averageBeat2,g1data.outlier,str2num(get(handles.cmin2,'String')),str2num(get(handles.cmax2,'String')),get(g1data.tfilt,'Value'),str2num(get(handles.beforeGUI,'String')),get(handles.apdbl,'Value'),str2num(get(handles.apdblnum,'String')),get(g1data.medifilt,'Value'));
%% Dual Maps 
if get(handles.diffopt,'Value') == 1
handles.acttime=handles.actmap1-handles.actmap2;
handles.decaytime=handles.mapR1-handles.mapR2;
handles.duration=handles.apmap1-handles.apmap2;
end

if get(handles.diffopt,'Value') == 2
handles.acttime=handles.actmap2-handles.actmap1;
handles.decaytime=handles.mapR2-handles.mapR1;
handles.duration=handles.apmap2-handles.apmap1;
end
%% Get rid of 0-0 values in dual maps
[rows cols]=size(handles.actmap1)
for r=1:rows
    for c=1:cols
        if handles.actmap1(r,c)==0 || handles.actmap2(r,c)==0
          handles.acttime(r,c)=NaN;
        end
                if handles.mapR1(r,c)==0 || handles.mapR2(r,c)==0
          handles.decaytime(r,c)=NaN;
                end
                if handles.apmap1(r,c)==0 || handles.apmap2(r,c)==0
          handles.duration(r,c)=NaN;
                end 
    end
end

% Update handles structure
act_times=handles.acttime(isnan(handles.acttime)==0);
dec_times=handles.decaytime(isnan(handles.decaytime)==0);
dur_times=handles.duration(isnan(handles.duration)==0);

pact=prctile(act_times,[5,50,95]);
pdec=prctile(dec_times,[5,50,95]);
pdur=prctile(dur_times,[5,50,95]);


% Results table
handles.rdata=zeros(3,4);
handles.rdata(1,1)=mean(act_times);handles.rdata(1,2)=std(act_times);handles.rdata(1,3)=std(act_times)/numel(act_times);handles.rdata(1,4)=((pact(3)-pact(1))/pact(2));
handles.rdata(2,1)=mean(dec_times);handles.rdata(2,2)=std(dec_times);handles.rdata(2,3)=std(dec_times)/numel(dec_times);handles.rdata(2,4)=((pdec(3)-pdec(1))/pdec(2));
handles.rdata(3,1)=mean(dur_times);handles.rdata(3,2)=std(dur_times);handles.rdata(3,3)=std(dur_times)/numel(dur_times);handles.rdata(3,4)=((pdur(3)-pdur(1))/pdur(2));
set(handles.rtable,'Data',handles.rdata);


guidata(hObject, handles);
zoom off
% Update handles structure

guidata(hObject, handles);
imchoice_Callback(hObject, eventdata, handles)
dualchoice_Callback(hObject, eventdata, handles)
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


% --- Executes on selection change in imchoice.
function imchoice_Callback(hObject, eventdata, handles)
% hObject    handle to imchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = findobj('Tag','ElectroMap');
g1data = guidata(h);
imchoice=get(handles.imchoice, 'Value');
%% Raw Images
if imchoice == 1
    axes(handles.map1)
    cla
    imshow(handles.frame1,[],'InitialMagnification', 400)
    colormap('gray')
    freezeColors
    hold on
    for i=1:size(handles.boundaries,1)
        plot(handles.boundaries{i}(:,2),handles.boundaries{i}(:,1),'b','LineWidth',2);
    end
    hold off
    axes(handles.cb1);
    cla reset
    hcb=colorbar;
    hcb.Location = 'southoutside'
    ax = gca;
    axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    axis off
    axes(handles.map2)
    cla
    imshow(handles.frame12,[],'InitialMagnification', 400)
    colormap('gray')
    freezeColors
    hold on
    for i=1:size(handles.boundaries,1)
        plot(handles.boundaries{i}(:,2),handles.boundaries{i}(:,1),'r','LineWidth',2);
    end
    hold off
        axes(handles.cb2);
    cla reset
    hcb=colorbar;
    hcb.Location = 'southoutside'
    ax = gca;
    %axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    axis off
end

%% Activation Maps
if imchoice == 2
axes(handles.map1)
cla
if get(handles.apdscale1,'Value') == 1
mini=min(min(handles.actmap1(handles.actmap1>0)));
maxi=max(max(handles.actmap1(handles.actmap1>0)));
end
if get(handles.apdscale1,'Value') == 2
mini=str2num(get(handles.cmin1,'String'));
maxi=str2num(get(handles.cmax1,'String'));
end
imshow(handles.actmap1, [mini maxi], 'InitialMagnification', 800),
    hold on
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    caxis([mini maxi]);
     freezeColors

    % colorbar
    axes(handles.cb1);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    %axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.apdscale1,'Value') == 1
        stepp=(maxi-mini)/5;
        hcb.TickLabels=[floor(mini):stepp:ceil(maxi)];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Activation Time (ms)';
        axis off
    end
    if get(handles.apdscale1,'Value') == 2
        stepp=(str2num(get(handles.cmax1,'String'))-str2num(get(handles.cmin1,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin1,'String')):stepp:str2num(get(handles.cmax1,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Activation Time (ms)';
        axis off
    end

axes(handles.map2)
cla
if get(handles.apdscale2,'Value') == 1
mini=min(min(handles.actmap2(handles.actmap2>0)));
maxi=max(max(handles.actmap2(handles.actmap2>0)));
end
if get(handles.apdscale2,'Value') == 2
mini=str2num(get(handles.cmin2,'String'));
maxi=str2num(get(handles.cmax2,'String'));
end
imshow(handles.actmap2, [mini maxi], 'InitialMagnification', 800),
    hold on
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    caxis([mini maxi]);
    freezeColors

% colorbar
    axes(handles.cb2);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    %axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.apdscale2,'Value') == 1
        stepp=(maxi-mini)/5;
        hcb.TickLabels=[floor(mini):stepp:ceil(maxi)];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Activation Time (ms)';
        axis off
    end
    if get(handles.apdscale2,'Value') == 2
        stepp=(str2num(get(handles.cmax2,'String'))-str2num(get(handles.cmin2,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin2,'String')):stepp:str2num(get(handles.cmax2,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Activation Time (ms)';
        axis off
end
end

%% APD/CaD Maps
if imchoice == 3
   axes(handles.map1)
   cla
    if get(handles.apdscale1,'Value') == 1
apmini1=min(min(handles.apmap1(handles.apmap1>0)));
apmaxi1=max(max(handles.apmap1(handles.apmap1>0)));
end
if get(handles.apdscale1,'Value') == 2
apmini1=str2num(get(handles.cmin1,'String'));
apmaxi1=str2num(get(handles.cmax1,'String'));
end
    imshow(handles.apmap1,[apmini1 apmaxi1], 'InitialMagnification', 800);
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    freezeColors
    axes(handles.cb1);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    %axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.apdscale1,'Value') == 1
        stepp=(apmaxi1-apmini1)/5;
        hcb.TickLabels=[floor(apmini1):stepp:ceil(apmaxi1)];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Duaration (ms)';
        axis off
    end
    if get(handles.apdscale1,'Value') == 2
        stepp=(str2num(get(handles.cmax1,'String'))-str2num(get(handles.cmin1,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin1,'String')):stepp:str2num(get(handles.cmax1,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Duartion (ms)';
        axis off
end       


    axes(handles.map2)
    cla
        if get(handles.apdscale2,'Value') == 1
apmini2=min(min(handles.apmap2(handles.actmap2>0)));
apmaxi2=max(max(handles.apmap2(handles.actmap2>0)));
end
if get(handles.apdscale2,'Value') == 2
apmini2=str2num(get(handles.cmin2,'String'));
apmaxi2=str2num(get(handles.cmax2,'String'));
end
    imshow(handles.apmap2,[apmini2 apmaxi2], 'InitialMagnification', 800);
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    freezeColors
    axes(handles.cb2);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    %axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.apdscale2,'Value') == 1
        stepp=(apmaxi2-apmini2)/5;
        hcb.TickLabels=[floor(apmini2):stepp:ceil(apmaxi2)];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Duaration (ms)';
        axis off
    end
    if get(handles.apdscale2,'Value') == 2
        stepp=(str2num(get(handles.cmax2,'String'))-str2num(get(handles.cmin2,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmin2,'String')):stepp:str2num(get(handles.cmax2,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='Duartion (ms)';
        axis off
end  
end
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns imchoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imchoice


% --- Executes during object creation, after setting all properties.
function imchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dualchoice.
function dualchoice_Callback(hObject, eventdata, handles)
% hObject    handle to dualchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
dualchoice=get(handles.dualchoice, 'Value');
if dualchoice == 1
    axes(handles.map3)
    cla
if get(handles.dualscale,'Value') == 1
    actmini=min(min(handles.acttime(isnan(handles.acttime)==0)));
    actmaxi=max(max(handles.acttime(isnan(handles.acttime)==0)));
end
if get(handles.dualscale,'Value') == 2
actmini=str2num(get(handles.cmind,'String'));
actmaxi=str2num(get(handles.cmaxd,'String'));
end
    imshow(handles.acttime,[actmini actmaxi], 'InitialMagnification', 800);
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    %jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    freezeColors
    axes(handles.cbd);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    %axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.dualscale,'Value') == 1
        stepp=(actmaxi-actmini)/5;
        hcb.TickLabels=[floor(actmini):stepp:ceil(actmaxi)];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='\Delta Activation Time (ms)';
        axis off
    end
    if get(handles.dualscale,'Value') == 2
        stepp=(str2num(get(handles.cmaxd,'String'))-str2num(get(handles.cmind,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmind,'String')):stepp:str2num(get(handles.cmaxd,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='\Delta Activation Time (ms)';
        axis off
        
end 
    
end

if dualchoice == 2
    axes(handles.map3)
    cla
    if get(handles.dualscale,'Value') == 1
    decmini=min(min(handles.decaytime(isnan(handles.decaytime)==0)))
    decmaxi=max(max(handles.decaytime(isnan(handles.decaytime)==0)))
    end
    if get(handles.dualscale,'Value') == 2
    decmini=str2num(get(handles.cmind,'String'));
    decmaxi=str2num(get(handles.cmaxd,'String'));
     end

    if decmini == 0
       decmini=min(min(handles.decaytime(handles.decaytime>0)))
    end
    imshow(handles.decaytime,[decmini decmaxi], 'InitialMagnification', 800);
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    %jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    freezeColors
    freezeColors
    axes(handles.cbd);
    cla reset
    hcb=colorbar;
    pretty=get(handles.colmap,'String'); jetcolormap = (colormap(pretty{get(handles.colmap,'Value')}));
    hcb.Location = 'southoutside'
    ax = gca;
    %axpos = ax.Position;
    cpos = hcb.Position;
    cpos(4) = 4*cpos(4);
    hcb.Position = cpos;
    %ax.Position = axpos;
    hcb.TicksMode='manual';hcb.TickLabelsMode='manual';
    if get(handles.dualscale,'Value') == 1
        stepp=(decmaxi-decmini)/5;
        hcb.TickLabels=[floor(decmini):stepp:ceil(decmaxi)];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='\Delta Duartion (ms)';
        axis off
    end
    if get(handles.dualscale,'Value') == 2
        stepp=(str2num(get(handles.cmaxd,'String'))-str2num(get(handles.cmind,'String')))/5;
        hcb.TickLabels=[str2num(get(handles.cmind,'String')):stepp:str2num(get(handles.cmaxd,'String'))];
        hcb.Ticks=[0.01,0.2:0.2:1];
        hcb.Label.String='\Delta Duartion (ms)';
        axis off
end 
end

if dualchoice == 3
    axes(handles.map3)
    cla
     handles.duration
    durmini=min(min(handles.duration(isnan(handles.duration)==0)))
    durmaxi=max(max(handles.duration(isnan(handles.duration)==0)))
    if durmini == 0
       durmini=min(min(handles.duration(handles.duration>0)))
    end
    imshow(handles.duration,[durmini durmaxi], 'InitialMagnification', 800);
    pretty=get(g1data.colmap,'String'); jetcolormap = (colormap(pretty{get(g1data.colmap,'Value')}));
    %jetcolormap(1,:) = [1, 1, 1];
    colormap(jetcolormap);
    freezeColors
end
% Hints: contents = cellstr(get(hObject,'String')) returns dualchoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dualchoice


% --- Executes during object creation, after setting all properties.
function dualchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dualchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in process2.
function process2_Callback(hObject, eventdata, handles)
% hObject    handle to process2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in process1.
function process1_Callback(hObject, eventdata, handles)
% hObject    handle to process1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function getMousePositionOnImage(hObject, eventdata, handles)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
cursorPoint = get(handles.map3, 'CurrentPoint');
curX = cursorPoint(1,1);
curY = cursorPoint(1,2);
handles.row=floor(curY);
handles.col=floor(curX);
[rows,cols,num] = size(g1data.images(:,:,:));
xLimits = get(handles.map3, 'xlim');
yLimits = get(handles.map3, 'ylim');
axes(handles.map3)
if rows == 1
        row = 1 
        curY= min(yLimits)+((max(yLimits)-min(yLimits))/2);
    end
    if cols == 1
        col = 1 
        curX= min(xLimits)+((max(xLimits)-min(xLimits))/2);
    end 

if (curX > min(xLimits) && curX < max(xLimits) && curY > min(yLimits) && curY < max(yLimits))
    disp('Hi')
    if isempty(handles.row) == 0
cla
dualchoice_Callback(hObject, eventdata, handles)
axes(handles.map3)
   hold on
plot(handles.col,handles.row,'k+','MarkerSize',20,'LineWidth',5);
hold off
    end
    axes(handles.dualaxes)
    cla
   [handles.sig1]=mapsbabyonepix_dual(get(handles.velalgo,'Value'),str2num(get(g1data.framerate,'String')),str2num(get(g1data.t,'String')),g1data.I,handles.images1,handles.averageBeat1,handles.row,handles.col,'b',str2num(get(g1data.beforeGUI,'String')),str2num(get(g1data.afterGUI,'String')),get(g1data.apdbl,'Value'),str2num(get(g1data.apdblnum,'String')),str2num(get(g1data.taustart,'String')),str2num(get(g1data.taufinish,'String')),1)
   hold on
   [handles.sig2]=mapsbabyonepix_dual(get(handles.velalgo,'Value'),str2num(get(g1data.framerate,'String')),str2num(get(g1data.t,'String')),g1data.I,handles.images2,handles.averageBeat2,handles.row,handles.col,'r',str2num(get(g1data.beforeGUI,'String')),str2num(get(g1data.afterGUI,'String')),get(g1data.apdbl,'Value'),str2num(get(g1data.apdblnum,'String')),str2num(get(g1data.taustart,'String')),str2num(get(g1data.taufinish,'String')),1)
end
guidata(hObject, handles);

% --- Executes on selection change in BLopt1.
function BLopt1_Callback(hObject, eventdata, handles)
% hObject    handle to BLopt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BLopt1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BLopt1


% --- Executes during object creation, after setting all properties.
function BLopt1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BLopt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thlen1_Callback(hObject, eventdata, handles)
% hObject    handle to thlen1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thlen1 as text
%        str2double(get(hObject,'String')) returns contents of thlen1 as a double


% --- Executes during object creation, after setting all properties.
function thlen1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thlen1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tfilt1.
function tfilt1_Callback(hObject, eventdata, handles)
% hObject    handle to tfilt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tfilt1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tfilt1


% --- Executes during object creation, after setting all properties.
function tfilt1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfilt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in usespline1.
function usespline1_Callback(hObject, eventdata, handles)
% hObject    handle to usespline1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns usespline1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from usespline1


% --- Executes during object creation, after setting all properties.
function usespline1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usespline1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function splineN1_Callback(hObject, eventdata, handles)
% hObject    handle to splineN1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of splineN1 as text
%        str2double(get(hObject,'String')) returns contents of splineN1 as a double


% --- Executes during object creation, after setting all properties.
function splineN1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to splineN1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function apdblnum1_Callback(hObject, eventdata, handles)
% hObject    handle to apdblnum1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of apdblnum1 as text
%        str2double(get(hObject,'String')) returns contents of apdblnum1 as a double


% --- Executes during object creation, after setting all properties.
function apdblnum1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdblnum1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apdbl1.
function apdbl1_Callback(hObject, eventdata, handles)
% hObject    handle to apdbl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns apdbl1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apdbl1


% --- Executes during object creation, after setting all properties.
function apdbl1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdbl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes on selection change in sfilt1.
function sfilt1_Callback(hObject, eventdata, handles)
% hObject    handle to sfilt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sfilt1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sfilt1


% --- Executes during object creation, after setting all properties.
function sfilt1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfilt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sfiltsize1_Callback(hObject, eventdata, handles)
% hObject    handle to sfiltsize1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sfiltsize1 as text
%        str2double(get(hObject,'String')) returns contents of sfiltsize1 as a double


% --- Executes during object creation, after setting all properties.
function sfiltsize1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfiltsize1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sfiltsigma1_Callback(hObject, eventdata, handles)
% hObject    handle to sfiltsigma1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sfiltsigma1 as text
%        str2double(get(hObject,'String')) returns contents of sfiltsigma1 as a double


% --- Executes during object creation, after setting all properties.
function sfiltsigma1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfiltsigma1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in velalgo1.
function velalgo1_Callback(hObject, eventdata, handles)
% hObject    handle to velalgo1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns velalgo1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velalgo1


% --- Executes during object creation, after setting all properties.
function velalgo1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velalgo1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t1_Callback(hObject, eventdata, handles)
% hObject    handle to t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t1 as text
%        str2double(get(hObject,'String')) returns contents of t1 as a double


% --- Executes during object creation, after setting all properties.
function t1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in invertopt1.
function invertopt1_Callback(hObject, eventdata, handles)
% hObject    handle to invertopt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invertopt1


% --- Executes on button press in invertopt2.
function invertopt2_Callback(hObject, eventdata, handles)
% hObject    handle to invertopt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invertopt2



function t2_Callback(hObject, eventdata, handles)
% hObject    handle to t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2 as text
%        str2double(get(hObject,'String')) returns contents of t2 as a double


% --- Executes during object creation, after setting all properties.
function t2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in velalgo2.
function velalgo2_Callback(hObject, eventdata, handles)
% hObject    handle to velalgo2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns velalgo2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velalgo2


% --- Executes during object creation, after setting all properties.
function velalgo2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velalgo2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sfiltsigma2_Callback(hObject, eventdata, handles)
% hObject    handle to sfiltsigma2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sfiltsigma2 as text
%        str2double(get(hObject,'String')) returns contents of sfiltsigma2 as a double


% --- Executes during object creation, after setting all properties.
function sfiltsigma2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfiltsigma2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sfiltsize2_Callback(hObject, eventdata, handles)
% hObject    handle to sfiltsize2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sfiltsize2 as text
%        str2double(get(hObject,'String')) returns contents of sfiltsize2 as a double


% --- Executes during object creation, after setting all properties.
function sfiltsize2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfiltsize2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sfilt2.
function sfilt2_Callback(hObject, eventdata, handles)
% hObject    handle to sfilt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sfilt2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sfilt2


% --- Executes during object creation, after setting all properties.
function sfilt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfilt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in aptime2.
function aptime2_Callback(hObject, eventdata, handles)
% hObject    handle to aptime2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns aptime2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from aptime2


% --- Executes during object creation, after setting all properties.
function aptime2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aptime2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apdbl2.
function apdbl2_Callback(hObject, eventdata, handles)
% hObject    handle to apdbl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns apdbl2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apdbl2


% --- Executes during object creation, after setting all properties.
function apdbl2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdbl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function apdblnum2_Callback(hObject, eventdata, handles)
% hObject    handle to apdblnum2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of apdblnum2 as text
%        str2double(get(hObject,'String')) returns contents of apdblnum2 as a double


% --- Executes during object creation, after setting all properties.
function apdblnum2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdblnum2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function splineN2_Callback(hObject, eventdata, handles)
% hObject    handle to splineN2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of splineN2 as text
%        str2double(get(hObject,'String')) returns contents of splineN2 as a double


% --- Executes during object creation, after setting all properties.
function splineN2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to splineN2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in usespline2.
function usespline2_Callback(hObject, eventdata, handles)
% hObject    handle to usespline2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns usespline2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from usespline2


% --- Executes during object creation, after setting all properties.
function usespline2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usespline2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tfilt2.
function tfilt2_Callback(hObject, eventdata, handles)
% hObject    handle to tfilt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tfilt2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tfilt2


% --- Executes during object creation, after setting all properties.
function tfilt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfilt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thlen2_Callback(hObject, eventdata, handles)
% hObject    handle to thlen2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thlen2 as text
%        str2double(get(hObject,'String')) returns contents of thlen2 as a double


% --- Executes during object creation, after setting all properties.
function thlen2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thlen2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BLopt2.
function BLopt2_Callback(hObject, eventdata, handles)
% hObject    handle to BLopt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BLopt2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BLopt2


% --- Executes during object creation, after setting all properties.
function BLopt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BLopt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in process2.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to process2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in zoomon.
function zoomon_Callback(hObject, eventdata, handles)
% hObject    handle to zoomon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoomon



function cmin1_Callback(hObject, eventdata, handles)
% hObject    handle to cmin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmin1 as text
%        str2double(get(hObject,'String')) returns contents of cmin1 as a double


% --- Executes during object creation, after setting all properties.
function cmin1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmax1_Callback(hObject, eventdata, handles)
% hObject    handle to cmax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmax1 as text
%        str2double(get(hObject,'String')) returns contents of cmax1 as a double


% --- Executes during object creation, after setting all properties.
function cmax1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apdscale1.
function apdscale1_Callback(hObject, eventdata, handles)
% hObject    handle to apdscale1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns apdscale1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apdscale1


% --- Executes during object creation, after setting all properties.
function apdscale1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdscale1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmin2_Callback(hObject, eventdata, handles)
% hObject    handle to cmin2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmin2 as text
%        str2double(get(hObject,'String')) returns contents of cmin2 as a double


% --- Executes during object creation, after setting all properties.
function cmin2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmin2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmax2_Callback(hObject, eventdata, handles)
% hObject    handle to cmax2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmax2 as text
%        str2double(get(hObject,'String')) returns contents of cmax2 as a double


% --- Executes during object creation, after setting all properties.
function cmax2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmax2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apdscale2.
function apdscale2_Callback(hObject, eventdata, handles)
% hObject    handle to apdscale2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns apdscale2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apdscale2


% --- Executes during object creation, after setting all properties.
function apdscale2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apdscale2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmind_Callback(hObject, eventdata, handles)
% hObject    handle to cmind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmind as text
%        str2double(get(hObject,'String')) returns contents of cmind as a double


% --- Executes during object creation, after setting all properties.
function cmind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmaxd_Callback(hObject, eventdata, handles)
% hObject    handle to cmaxd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmaxd as text
%        str2double(get(hObject,'String')) returns contents of cmaxd as a double


% --- Executes during object creation, after setting all properties.
function cmaxd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmaxd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dualscale.
function dualscale_Callback(hObject, eventdata, handles)
% hObject    handle to dualscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dualscale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dualscale


% --- Executes during object creation, after setting all properties.
function dualscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dualscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in drawcon1.
function drawcon1_Callback(hObject, eventdata, handles)
% hObject    handle to drawcon1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of drawcon1



function conbound1_Callback(hObject, eventdata, handles)
% hObject    handle to conbound1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of conbound1 as text
%        str2double(get(hObject,'String')) returns contents of conbound1 as a double


% --- Executes during object creation, after setting all properties.
function conbound1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conbound1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in drawcon2.
function drawcon2_Callback(hObject, eventdata, handles)
% hObject    handle to drawcon2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of drawcon2



function conbound2_Callback(hObject, eventdata, handles)
% hObject    handle to conbound2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of conbound2 as text
%        str2double(get(hObject,'String')) returns contents of conbound2 as a double


% --- Executes during object creation, after setting all properties.
function conbound2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conbound2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in drawcond.
function drawcond_Callback(hObject, eventdata, handles)
% hObject    handle to drawcond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of drawcond



function conboundd_Callback(hObject, eventdata, handles)
% hObject    handle to conboundd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of conboundd as text
%        str2double(get(hObject,'String')) returns contents of conboundd as a double


% --- Executes during object creation, after setting all properties.
function conboundd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conboundd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in diffopt.
function diffopt_Callback(hObject, eventdata, handles)
% hObject    handle to diffopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns diffopt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from diffopt


% --- Executes during object creation, after setting all properties.
function diffopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diffopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mapapply.
function mapapply_Callback(hObject, eventdata, handles)
% hObject    handle to mapapply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

listbox2_Callback(hObject, eventdata, handles)
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in colmap.
function colmap_Callback(hObject, eventdata, handles)
% hObject    handle to colmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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


% --- Executes on selection change in velalgo.
function velalgo_Callback(hObject, eventdata, handles)
% hObject    handle to velalgo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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



function dura_Callback(hObject, eventdata, handles)
% hObject    handle to dura (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dura as text
%        str2double(get(hObject,'String')) returns contents of dura as a double


% --- Executes during object creation, after setting all properties.
function dura_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dura (see GCBO)
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



function dura2_Callback(hObject, eventdata, handles)
% hObject    handle to dura2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dura2 as text
%        str2double(get(hObject,'String')) returns contents of dura2 as a double


% --- Executes during object creation, after setting all properties.
function dura2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dura2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes during object deletion, before destroying properties.
function rtable_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to rtable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in sigexport.
function sigexport_Callback(hObject, eventdata, handles)
% hObject    handle to sigexport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
h = findobj('Tag','ElectroMap');
g1data = guidata(h);
[filename,pathname] = uiputfile({'*.csv';'*.txt';'*.mat'}, 'Save signal');
[~,~,ext] = fileparts(filename);
file=[pathname,filename];
handles
exposure=1/str2num(get(g1data.framerate,'String'));
tictoc=0:exposure:(length(handles.sig1)-1)*exposure;
T=table(tictoc',handles.sig1,handles.sig2);
writetable(T,file,'Delimiter',',','WriteVariableNames',false);