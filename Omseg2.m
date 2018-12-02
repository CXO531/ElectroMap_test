function [locs,minimas,q2locs,avgCL,numofpeaksoverall,newpeakheight]=Omseg2(signal,peakheight,peakdist,minpeakheight,minpeakdist,minboundary,segchoice,minmumofpeaks,numimages,div)
%% Function for segmenting a signal based on CL according to user settings
% Chris O'Shea 2018

%% Detect Peaks and minimums
normsig=signal-min(signal);
normsig=normsig./max(normsig);
insig=imcomplement(normsig);
insig=insig-min(insig);
insig=smooth(insig);
normsig
%find peaks
maxfluo=max(signal);
newpeakheight=maxfluo*peakheight;
[pks, locs] = findpeaks(signal, 'MINPEAKHEIGHT',newpeakheight, 'MINPEAKDISTANCE', peakdist);
remove1=0;
removeend=0;

%% find minimas
alocs=[];blocs=[];
bisg=[];asig=[];
minimas=zeros(length(locs),2);


if isempty (locs) == 0
bsig=insig(1:locs(1));
else
    locs=[1 2 3 4 5 6 7 8];
    bsig=[1 2 3 4 5 6];
end
if numel(bsig) < 3
    bsig = [1 2 3 4 5 6 7];
end
bsig
[ipks, blocs] = findpeaks(bsig, 'MINPEAKHEIGHT',0.5*peakheight*max(insig), 'MINPEAKDISTANCE', 1);
if isempty(blocs) ~= 1 
minimas(1,1)=blocs(end);
if numel(locs) > 1
asig=insig(locs(1):locs(2));
else
    asig=[1 2 3 4 5 6 7];
end
[ipks, alocs] = findpeaks(asig, 'MINPEAKHEIGHT',0.5*peakheight*max(insig), 'MINPEAKDISTANCE', 1);
if isempty(alocs) == 1
   [ipks,alocs]=max(asig);
end
minimas(1,2)=alocs(1)+locs(1);
loopstart=2
elseif isempty(blocs) == 1  %if no minimum before first peak, then ignore
    remove1=1;
end

for i=2:length(locs)-1
alocs=[];blocs=[];
bsig=insig(locs(i-1):locs(i));
if length(bsig) > 2
[ipks, blocs] = findpeaks(bsig, 'MINPEAKHEIGHT',0.5*peakheight*max(insig), 'MINPEAKDISTANCE', 2);
end
if isempty(blocs) == 1
    isempty(asig) == 0
    [ipks, blocs]=max(asig)
else
    asig=[1:20]
    [ipks, blocs]=max(asig)
end
if isempty(blocs) == 1
    blocs = 10
end
minimas(i,1)=blocs(end)+locs(i-1);
asig=insig(locs(i):locs(i+1));
if length(asig) > 2
[ipks, alocs] = findpeaks(asig, 'MINPEAKHEIGHT',0.5*peakheight*max(insig), 'MINPEAKDISTANCE', 2);
end
if isempty(alocs) == 1
    [ipks, alocs]=max(asig)
end
minimas(i,2)=alocs(1)+locs(i);
end    


alocs=[];blocs=[];
if numel(locs) > 1
bsig=insig(locs(numel(locs)-1):locs(numel(locs)));
else
    bsig=[1 2 3 4 5 6 7]
end
if numel(asig) > 2
[ipks, blocs] = findpeaks(bsig, 'MINPEAKHEIGHT',0.001*max(insig), 'MINPEAKDISTANCE', 2);
end
asig=insig(locs(numel(locs)):end);
if numel(asig) > 3
[ipks, alocs] = findpeaks(asig, 'MINPEAKHEIGHT',0.001*max(insig), 'MINPEAKDISTANCE', 2);
end
if isempty(blocs) ~= 1 && isempty(alocs) ~= 1
minimas(numel(locs),1)=blocs(1)+locs(numel(locs)-1);
minimas(numel(locs),2)=alocs(1)+locs(numel(locs));
else
    removeend=1
end

if remove1==1
   locs=locs(2:end)
   minimas=minimas(2:end,:) 
end

if removeend==1
   locs=locs(1:end-1)
   minimas=minimas(1:end-1,:) 
end

%% Split up 'averages'
CL=zeros((length(locs)),3);
CL(:,1)=locs(:);
CLdiff=zeros(1,(length(CL)-1));
%Calculate CLs
if length(locs) >= 2
CL(1,2)=locs(2)-locs(1) %first peak as no refrence before, so use next peak
end
for i=2:(length(locs))
    CL(i,2)=locs(i)-locs(i-1);
end

for i=2:(length(locs)-2)
    CLdiff(i)=(locs(i+1)-locs(i))-(locs(i)-locs(i-1));
    if abs(CLdiff(i)) > minboundary
        CL(i,3)=1;
        if CL(i-1,3) == 1
            CL(i,3) = 0
        end
    end
end

%Find diffrenet regions of constant CL , and position in signal.
q1=[];k=1;j=0; q2=[]; q2locs=[];
if isempty(locs) == 1 %no peaks
    h=errordlg('No action potential detected')
    waitfor(h)
end

if length(locs) > 1
    for i=1:(length(locs))
        if CL(i,3) == 0
            q1(k,i)=[CL(i,2)];
        else q1(k,i)=[CL(i,2)];
            k=k+1;
        end
    end
     
    for k=1:length(q1(:,2))
        numofpeaks=sum(q1(k,:)~=0);
        if numofpeaks >= minmumofpeaks
           j=j+1;
           q2(j,:)=q1(k,:);
        end
        numofpeaksoverall = numofpeaks;
    end
    
    if isempty(q2) ==1
        errordlg('Did not find enough peaks, please reduce minimum number required in Signal Processing Options. If signal possible arrhythmia suggest min number and min boundary both reduced to 1');
    end
     avgCL=[];
    for j=1:length(q2(1,:))
        for i=1:length(q2(:,1))
            if q2(i,j) ~=0
                q2locs(i,j)=locs(j);
            else
                q2locs(i,j)=0;
            end
        end
    end
    avgCL=zeros(2,length(q2(:,1)));
    for i=1:length(q2(:,1))
        pos=mean(find(q2(i,:)));
        %avgCL(1,i)=round(pos);
        avgCL(1,i)=(pos);
        avgCL(2,i) = mean(nonzeros(q2(i,:)));
    end
end    

%% check if first or last peak too close to end
% if q2locs(1,1) < bframes
%    q2locs=q2locs(:,2:end);
% end
% if q2locs(length(q2locs(:,1)),length(q2locs(1,:)))+aframes > num_images
%     q2locs=q2locs(:,1:end-1);
% end

%% Sub-Segmentation
if segchoice == 2 || segchoice == 3
    q2locsnew=[];
    k=0;col=0;row=1;
    for i=1:length(q2locs(:,1))
        for j=1:length(q2locs(1,:))
            if q2locs(i,j) ~= 0;
                k=k+1;
                colum=mod(k,div);
                
                if colum == 0
                    col=div;
                else col=colum;
                end
                
                q2locsnew(row,col)=q2locs(i,j);
                if j == length(q2(1,:)) %if loop to stop martix dim error
                    j=length(q2(1,:))-1;
                end
                
                
                if col == div
                    row=row+1;
                    k=0;
                end
                if j+1 < length(q2locs(1,:))
                    if q2locs(i,j+1)==0
                        row = row+2;
                    end
                end
            end
        end
    end
    
    
    if div > size(q2locsnew,2)
        div=size(q2locsnew,2)
    end
        
    q2locsnew=[q2locsnew;zeros(2,div)]; %added to allow 'buffer' at end to allow idneftication of last constat CL region.
    firstrow=1; A=[]; B=[];
    
    for i=1:length(q2locsnew(:,1))
        check=any((q2locsnew(i,:)));
        if check == 0
            A=q2locsnew(firstrow:i,:);
            A=sort(A,1);
            A=sort(A,2);
            B=[B;A];
            firstrow=i+1;
        end
    end
    
    q2locs=B(any(B,2),:);
    
     if segchoice == 3
        row =1; C=[];
        for i = 1:(length(B(:,1))-1)
            if any(B(i+1,:)) == 0
                C(row,:)=B(i,:);
                row = row+1;
            end
        end
        q2locs= C(any(C,2),:);
     end
     
     %update avgcl for subseq
     for i=1:length(q2locs(:,1));
        if any((q2locs(i,:))) == 1
            pos=mean(find(q2locs(i,:)));
            %avgCL(1,i)=round(pos);
            avgCL(1,i)=pos;
            for j=1:(length(q2locs(1,:))-1) %%does nothing when 1 peak?
                disp('f')
                cld(i,j)=q2locs(i,j+1)-q2locs(i,j);
                if q2locs(i,j+1) == 0 || q2locs(i,j) == 0
                    cld(i,j)=0;
                end
            end
            
            if length(q2locs(1,:)) ~= 1
                avgCL(2,i) = mean(nonzeros(cld(i,:)));
            else
                if i < length(q2locs)
                    avgCL(2,i)= locs(i+1)-locs(i);
                else avgCL(2,i)=avgCL(2,i-1);
                end %to deal with last peak
            end
        end
    end
    if length((q2locs(1,:))) == 1 %%only one peak
        q=1;
        for i=1:length(CL(:,1));
            if q < length(q2locs)
                if CL(i,1) == q2locs(q); %%check this peak in settings
                    avgCL(2,q) = CL(i,2);
                    q=q+1;
                end
            end
        end
    end
    
    avgCL=avgCL(:,any(avgCL,1));  %columns
    %avloc=locs(avgCL(1,:))
end
