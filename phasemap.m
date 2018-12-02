 function phasemap(images, mask, refsig, framerate, minpeakdist, before, after)
%fuction for taking a processed image stack and making a phase map
tic
[rows, cols, num]= size(images(:,:,:));
phases=zeros(rows,cols,num);
refsig=refsig-min(refsig);
[pks locs] = findpeaks(refsig, 'MINPEAKHEIGHT', max(refsig)/2, 'MINPEAKDISTANCE', minpeakdist);
hold on
CL=locs(2)-locs(1);
B = 1/CL*ones(CL,1); %MOVING AVERAGE FILTER FOR SHIFT
before=round(before*framerate);
after=round(after*framerate);
beatcheck=zeros(1,length(B));
if locs(1) < before
    start = 2
else start = 1
end
for i=start:length(locs)
    for j=-before:after
        locs(i)+j
        beatcheck(locs(i)+j)=1;
    end
end
beatcheck
%% find phase maps
for row = 1:rows
    for col =1:cols
        for i=1:num
            pixelsignal(i)=images(row,col,i);
        end
        %pixelsignal=imcomplement(pixelsignal);
        out = filter(B,1,pixelsignal);    
        shiftedsignal=pixelsignal-out;
        shiftedpixelsignal=shiftedsignal(CL:length(shiftedsignal));
        hsig=hilbert(shiftedpixelsignal);
        for i =1:num-CL

            phases(row,col,i)=-1*angle(hsig(i));
           
            %phases(row,col,i)= imag(hsig(i));
        end
    end
end

figure,
for i =1:num
    A=phases(:,:,i);
    for row = 1:rows
    for col =1:cols
        if mask(row,col)==0
            A(row,col)=NaN;
        end
    end
    end
    ddd=[0 0 0;jet];
    imshow(A, [-pi pi],'ColorMap',ddd,'InitialMagnification', 800);
end
toc

%% find PSs (based on Tomii et al 2016)
V=zeros(rows,cols,num);
Vhold=zeros(rows,cols);
for frame=1:num;
    B=phasenorm(:,:,frame);
    varmatrix=zeros(9,9);
for row = 5:rows-5
    for col =5:cols-5
        Vij=0;
        %find local phases
        for krow=1:9
            for kcol=1:9
                theta=(B((row-5)+krow,(col-5)+kcol));
                x=exp(1i*theta);
                varmatrix(krow,kcol)=x;
                Vij=Vij+x;
            end
        end
    Vij=Vij/81;
    Vij=abs(Vij);
    Vhold(row,col)=1-Vij;
    V(:,:,frame)=Vhold;
    end
end
   [r,c] = find(Vhold==max(Vhold(:)));
   if max(Vhold(:)) > 0.7
       PSI(frame,1) = r;
       PSI(frame,2) = c;
   else
       PSI(frame,1)=0;
       PSI(frame,2)=0;
   end
end

plot(PSI(:,1),PSI(:,2))
figure,
for i=1:num
imshow(V(:,:,i), [0 1], 'Colormap',jet, 'InitialMagnification', 800),
end
        



