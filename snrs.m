
function [mapSNRr,mapSNRdb,allSNRr,allSNRdb]=SNRs(images,mask,before,after,tfilt);
order=3; framesize=11;
for r = 1:size(images,1)
    for c = 1:size(images,2)
        if mask(r,c) ~= 0
            signalav = imcomplement(squeeze((images(r,c,:))));
            signalav = (double(signalav));
            
            if tfilt == 2
                signalav = sgolayfilt(signalav, order,framesize);
            end
            signalav=signalav-min(signalav);
            [maxval, maxInd] = max(signalav);
            noiseinds=zeros(numel(signalav),1);
            
            for m=-before:after
                if (maxInd+m) > 1 && (maxInd+m) <= numel(signalav)
                    noiseinds(maxInd+m)=1;
                end
            end
            %noiseinds'
            %pause(10)
            for j=1:numel(signalav)
                if noiseinds(j) == 0
                    noise(j)=signalav(j);
                else
                    noise(j)=NaN;
                end
            end
%             noiseinds
%             pause(14)
            amp1(r,c)=(maxval-nanmean(noise));
noise1(r,c)=(nanstd(noise)-nanmean(noise));

% if r == 25 && c == 50
%     figure,
%     plot(signalav,'k')
%     hold on
%     plot(noise,'r')
%     plot(maxInd,maxval,'ob')
%     line([0 numel(signalav)], [maxval maxval])
%     line([0 numel(signalav)],[nanmean(noise)+nanstd(noise) nanmean(noise)+nanstd(noise)])
%     line([0 numel(signalav)],[nanmean(noise) nanmean(noise)])
%     noiseinds
%     pause(10)
%     figure,
% end
mapSNRr(r,c)=(maxval-nanmean(noise))/(nanstd(noise)); %changed because pks amp also shirtd by mean noise. 
mapSNRdb(r,c)=20*log10(mapSNRr(r,c));

% if r == 25 && c == 50
%     figure,
%     plot(signalav,'k')
%     hold on
%     plot(noise,'r')
%     plot(maxInd,maxval,'ob')
%     line([0 numel(signalav)], [maxval maxval])
%     line([0 numel(signalav)],[nanmean(noise)+nanstd(noise) nanmean(noise)+nanstd(noise)])
%     line([0 numel(signalav)],[nanmean(noise) nanmean(noise)])
%     title(['SNR = ',num2str(mapSNRr(r,c)),'(',num2str(mapSNRdb(r,c)),'db)'])
% end

        else
            noise1(r,c)=NaN;
            mapSNRr(r,c)=NaN;
            mapSNRdb(r,c)=NaN;
        end
    end
end

allSNRr=mapSNRr(isnan(mapSNRr)~=1);
allSNRdb=mapSNRdb(isnan(mapSNRdb)~=1);
