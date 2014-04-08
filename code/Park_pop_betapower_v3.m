function [BetaOut, BetaControl] = Park_pop_betapower_v3 (chData, basetime, infusetime)

    srate=1000; %downsampled sampling rate;
    timebinlength=srate/20; %if 10,=100ms per point; if 20, 50ms each points, if 50, 20ms each point
    betaPower=(chData.beta).^2;
    totaltime=length(chData.beta)*timebinlength/1000; %in seconds
    
    %%
    window=5; %in sec
    winstep=5; % in sec
    
    InfuseBeta=nan(1,floor((totaltime-window)/winstep+1));
    for i=1:floor((totaltime-window)/winstep+1) %%sweaping window of 100s long, and 10s apart.
        t=1+(i-1)*winstep*1000/timebinlength;
        InfuseBeta(i)=nanmedian(betaPower(t:t+window*1000/timebinlength-1));
    end
    
    basepoints=(basetime-window)/winstep+1;
    infusepoints=(infusetime-window)/winstep+1;
    totalpoints=length(InfuseBeta);
    
    Basebeta=nanmedian(InfuseBeta(1:basepoints));
    
    %%
    %%to find changes in betapower
    timetottest=50; %50sec to compare to baseline
    
    Betachangettest=nan(1,totalpoints-timetottest/winstep);
    Meanbetachange=nan(1,totalpoints-timetottest/winstep);
    incbeta=[];
    decbeta=[];
    for i=1:(totalpoints-timetottest/winstep)
        Betachangettest(i)=ttest(InfuseBeta(i:(i+timetottest/winstep-1)), Basebeta, 0.05, 'both'); %p=0.05 or 0.01; most cases is 0.01
        Meanbetachange(i)=nanmean(InfuseBeta(i:(i+timetottest/winstep-1)));
        if Meanbetachange(i)>Basebeta && Betachangettest(i)==1
            incbeta(end+1)=i;
        end
        
        if Meanbetachange(i)<Basebeta && Betachangettest(i)==1
            decbeta(end+1)=i;
        end
    end
    
    %%
    %% to eliminate false identified changes.
    for i=1:2
        if i==1
            tempbeta=incbeta;
        else
            tempbeta=decbeta;
        end
        
        Tempchange=zeros(1,length(Betachangettest));
        Tempchange(tempbeta)=1;
        
        stringent=[1 1 1 1 1 1 1 1]; %%5 consecutive windows, that is 8*winstep=40sec.
        
        TempBeta=conv(Tempchange, stringent, 'valid');
        TempBeta=horzcat(zeros(1, length(stringent)-1), TempBeta);
        
        TrueChange{i}=find(TempBeta>length(stringent)-2);
        if isempty (TrueChange{i})
            BetaOutPut(i,1)=NaN;  %latency of beta change onset
            BetaOutPut(i,2)=NaN;  %peak of beta power change
        else
            BetaOutPut(i,1)=TrueChange{i}(1)*winstep;
            if i==1
                BetaOutPut(i,2)=max(Meanbetachange(TrueChange{i}));
            end
            if i==2
                BetaOutPut(i,2)=min(Meanbetachange(TrueChange{i}));
            end
        end
        
    end
    
    recoverbeta=nanmedian(InfuseBeta(totalpoints-basepoints+1:totalpoints));
    
    BetaOut=[Basebeta BetaOutPut(1,:) BetaOutPut(2,:) recoverbeta];
    
    BetaControl=[Basebeta 750 nanmedian(InfuseBeta(150:150+basepoints-1)) NaN NaN recoverbeta];
    
    
    
    figure (54)
    plot(1:winstep:length(InfuseBeta)*winstep, InfuseBeta, 'r-')
    hold on;
    scatter(1:winstep:length(Betachangettest)*winstep, Betachangettest)
    
%     figure (1)
%     x=0:timebinlength/1000:length(chData.beta)*timebinlength/1000;
%     y=1:200;
%     imagesc(x, y, log(chData.allfreq))
%     axis xy
%     axis ([0 length(chData.beta)*timebinlength/1000 4 70])
%     xlabel ('time (s)');
%     ylabel ('Frequency (Hz)')
%     colorbar;
%     
end
