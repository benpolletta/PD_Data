function [] = parkpop(folderName, filetocombine, channels)

srate=1000; %downsampled sampling rate;
timebinlength=srate/20; %if 10,=100ms per point; if 20, 50ms each points, if 50, 20ms each point

% [DataDirectory, DataAnalysisDirectory] = SetDirectory_temp(folderName);

folderName = [folderName, '/'];

parfor f=1:length(filetocombine)

    data = load ([folderName, filetocombine{f}]);
    
    for a=1:channels
        
        channel=['ch' num2str(a)];
        totaltime=length(data.(channel).energy);
        downPoint=floor(totaltime/timebinlength);
        
        %energy=squeeze(energy);
        timepoint=size(data.(channel).energy,1);
        
        for i=1:downPoint
            boundry=(i-1)*timebinlength+1:i*timebinlength;
            
            tempallfreq=data.(channel).energy(1:timepoint, boundry)';
            tempphi=data.(channel).phase(1:timepoint,boundry)';
            tempbeta=data.(channel).betaenergy(boundry);
            tempLgamma=data.(channel).logammaenergy(boundry);
            tempHgamma=data.(channel).higammaenergy(boundry);
            tempLbeta=data.(channel).lobetaenergy(boundry);
            tempHbeta=data.(channel).hibetaenergy(boundry);
            
            popDataLocal(f).(channel).allfreq(1:timepoint,i)=  median(tempallfreq);
            popDataLocal(f).(channel).phi(1:timepoint,i)= median(tempphi);
            popDataLocal(f).(channel).beta(i)= median(tempbeta);
            popDataLocal(f).(channel).logamma(i)= median(tempLgamma);
            popDataLocal(f).(channel).higamma(i)= median(tempHgamma);
            popDataLocal(f).(channel).lobeta(i)= median(tempLbeta);
            popDataLocal(f).(channel).hibeta(i)=median(tempHbeta);
        end
        
        popDataLocal(f).(channel).LFP=data.(channel).LFP;
        %         clear energy boundry betaenergy logammaenergy higammaenergy totaltime downPoint i tempallfreq tempbeta tempHgamma tempLgamma lobetaenergy hibetaenergy tempHbeta tempLbeta;
    end
    
end

% fileName=[filetocombine(1,1:8),'_',filetocombine(end,6:8)];
outputname=[folderName,'ParkPop.mat'];

for a=1:channels
    chan=['ch' num2str(a)];
    tempData = popDataLocal(:).(chan);
    
    popData.(chan).allfreq=[];
    popData.(chan).phi=[];
    popData.(chan).beta=[];
    popData.(chan).logamma=[];
    popData.(chan).higamma=[];
    popData.(chan).lobeta=[];
    popData.(chan).hibeta=[];
    popData.(chan).LFP=[];
    
    popFields = fields(tempData);
    
    for n = 1:numel(popFields)
        pdLocalChan = cat(1,popDataLocal(:).(chan));
        popFieldName = popFields{n};
        popData.(chan).(popFieldName) = horzcat(pdLocalChan(:).(popFieldName)); % could also use cat(2,...
    end
    
end

popData.files=filetocombine;
save (outputname, '-struct','popData','-v6');

end
% %%
% load (outputname)
% 
% % end
% 
% figure (2)
% x=0:timebinlength/1000:length(Park.(fileName).beta)*timebinlength/1000;
% y=1:200;
% imagesc(x, y, log(Park.(fileName).allfreq))
% axis xy
% axis ([0 length(Park.(fileName).beta)*timebinlength/1000 4 200])
% xlabel ('time (s)');
% ylabel ('Frequency (Hz)')
% colorbar;
% saveas(gcf,[FigDirectory, fileName, '_ch',num2str(chan),'_spec.jpg'])
% saveas(gcf,[FigDirectory, fileName, '_ch',num2str(chan),'_spec.fig'])
% 
% figure (5)
% x=0:timebinlength/1000:length(Park.(fileName).beta)*timebinlength/1000;
% y=1:200;
% imagesc(x, y, Park.(fileName).phi)
% axis xy
% axis ([0 length(Park.(fileName).beta)*timebinlength/1000 4 70])
% xlabel ('time (s)');
% ylabel ('Phase')
% caxis([-pi pi])
% colorbar;
% saveas(gcf,[FigDirectory, fileName, '_ch',num2str(chan),'_phi.jpg'])
% saveas(gcf,[FigDirectory, fileName, '_ch',num2str(chan),'_phi.fig'])
% 
% figure (3)
% x=timebinlength/1000*(1:length(Park.(fileName).beta));
% line(x,log(Park.(fileName).beta))
% xlabel ('Time (s)');
% ylabel ('Log Power')
% title('Beta Power')
% saveas(gcf,[FigDirectory, fileName, '_ch',num2str(chan),'_power.jpg'])
% saveas(gcf,[FigDirectory, fileName, '_ch',num2str(chan),'_power.fig'])



% figure (4)
% x=(1:length(Park.(filename).LFP))/1000;
% plot(x,Park.(filename).LFP-smooth(Park.(filename).LFP,500)','k')
% xlabel ('Time (s)');
% ylabel ('Voltage (mV)')
% title('LFP')
%
% saveas(gcf,[FigDirectory, filename, '_ch',num2str(channel),'_filtLFP.jpg'])
% saveas(gcf,[FigDirectory, filename, '_ch',num2str(channel),'_filtLFP.fig'])

%%
% timeperiod=200*srate/timebinlength; %sec of data looked.
% beftime=50*srate/timebinlength;  %before infusion
% durtime=1200*srate/timebinlength; %time to start calculate during infusion
% recovertime=3000*srate/timebinlength;
% tempbeta=Park.(filename).beta;

%% then run this:


%betaPower=[nanmean(tempbeta(beftime:beftime+timeperiod)),nanmean(tempbeta(durtime:durtime+timeperiod)),nanmean(tempbeta(recovertime:recovertime+timeperiod))];

