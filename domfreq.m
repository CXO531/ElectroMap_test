function[freqmap] = domfreq(mask,imagestack,framerate,minf,maxf,fbin,winopt)
wb=waitbar(0.1,'Calculating Dominant Frequncies');
[rows,cols,num] = size(imagestack(:,:,:));
freqmap=zeros(rows,cols);
imagestack=double(imagestack);
framerate=framerate*1000;
T=(1/framerate);
Fs=framerate;
L=num;
t=(0:L-1)*T;
order = 3;
framesize = 11;
if winopt == 0
winvec=ones(num,1)
end
if winopt == 1
winvec=hann(num)';
end
% window = hann(T*framerate);
padf=(Fs)/(fbin*num)
lp = padf*(num);
lpa = nextpow2(lp);
lpad=2.^lpa;
%mask=reshape(mask,[numel(mask),1]);
%imagestack=reshape(imagestack,[rows*cols,num]);
for row = 1:rows
                row
    for col = 1:cols
        if mask(row,col)~=0
        pixelsignal=squeeze(imagestack(row,col,:))';
        %pixelsignal=imcomplement(pixelsignal);
        %pixelsignal = sgolayfilt(pixelsignal, order,framesize);
        pixelsignal = pixelsignal - mean(pixelsignal); 
        pixelsignal=pixelsignal.*winvec;
        pixelsignal=pixelsignal(1,:);
%         if row==10 && col==32
%           figure,
%           pixelsignal
%           plot(pixelsignal)
%         end    
%         
%                 if row==40 && col==32
%           figure,
%           pixelsignal
%           plot(pixelsignal)
%         end  
        Y=fft(pixelsignal,lpad);
        Y=Y(1:lpad/2+1);
        P2=abs(Y/L);
        P1 = P2;
        P1(2:end-1) = 2*P1(2:end-1);
%        if row==10 && col==32
%           figure,
%           plot(P1)
%        end 
%         
% %               if row==40 && col==32
% %           figure,
% %           plot(P1)
% %         end 
        %P1(1:20)=0;
        f = 0:(Fs/lpad):(Fs/2);
        [o,oi]=find(f<maxf);
        [u,ui]=find(f>minf);
        [maxpower,maxind]=max(P1(min(ui):max(oi)));
        [~,xx]=find(P1==maxpower);
        xx=xx(1);
        freqmap(row,col)=f(xx);
        else
           freqmap(row,col)=NaN;
        end
    end
end
delete(wb)
% for r=1:(rows*cols)
%         if mask(r)~=0
%         pixelsignal=squeeze(imagestack(r,:))';
%         %pixelsignal=imcomplement(pixelsignal);
%         %pixelsignal = sgolayfilt(pixelsignal, order,framesize);
%         pixelsignal = pixelsignal - mean(pixelsignal); 
%         pixelsignal=pixelsignal.*winvec;
%         Y=fft(pixelsignal,lpad);
%         Y=Y(1:lpad/2+1);
%         P2=abs(Y/L);
%         P1 = P2;
%         P1(2:end-1) = 2*P1(2:end-1);
%         %P1(1:20)=0;
%         f = 0:(Fs/lpad):(Fs/2);
%         [o,oi]=find(f<maxf);
%         [u,ui]=find(f>minf);
%         [maxpower,maxind]=max(P1(min(ui):max(oi)));
% %         [~,xx]=find(P1==maxpower)
% %         maxpower
% %         P1
%         [xx]=find(P1==maxpower)
%         freqmap(r)=f(xx);
%         else
%            freqmap(r)=NaN;
%         end
%     end
% freqmap=reshape(freqmap,[rows cols]);
% mask=reshape(mask,[rows cols]);

freqmap=double(freqmap);
mask=double(mask);
freqmap=freqmap.*mask;


% 
%         row=15;
%         col=15;
%         pixelsignal=squeeze(imagestack(row,col,:))';
%         %pixelsignal=imcomplement(pixelsignal);
%         %pixelsignal = sgolayfilt(pixelsignal, order,framesize);
%         pixelsignal = pixelsignal - mean(pixelsignal); 
%         figure,
%         plot(pixelsignal)
% %         hold on
% %         plot(pixelsignal.*winvec);
%         %winvec = hamming(length(pixelsignal));
%         Y=fft(pixelsignal,lpad);
%         Y=Y(1:lpad/2+1);
%         P2=abs(Y/L);
%         P1 = P2;
%         P1(2:end-1) = 2*P1(2:end-1);
%         %P1(1:20)=0;
%         f = 0:(Fs/lpad):(Fs/2);
%         [o,oi]=find(f<maxf);
%         [u,ui]=find(f>minf);
%         [maxpower,maxind]=max(P1(min(ui):max(oi)));
%         [~,xx]=find(P1==maxpower);
%         freqmap(row,col)=f(xx);
%         figure,
%         plot(f,P1)
%         figure,
