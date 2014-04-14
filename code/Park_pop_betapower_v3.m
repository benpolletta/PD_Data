function [BetaOut, BetaControl] = Park_pop_betapower_v3 (chData, basetime, infusetime)

    srate = 1000; % Downsampled sampling rate.
    binrate = 20; % How many timebins per second.
    timebinlength = srate/binrate; % If 10, 100ms per point; if 20, 50ms each points; if 50, 20ms each point.
    betaPower = (chData.beta).^2; % Squared magnitude of FFT; should be sampled at srate.
    totaltime = length(chData.beta)/binrate; % In seconds?
    
    %%
    window = 5; % In sec.
    winstep = 5; % In sec.
    
    winstep_data_pts = winstep*1000/timebinlength;
    window_data_pts = window*1000/timebinlength;
    no_windows = floor((totaltime-window)/winstep+1); % Total # windows, based on size of steps b/ween windows.
    
    InfuseBeta=nan(1,no_windows); % Vector to contain beta, length of # windows.
    for i=1:no_windows % Sweeping window of 100s long, and 10s apart.
        t=1+(i-1)*winstep_data_pts; % Window start time. 
        InfuseBeta(i)=nanmedian(betaPower(t:t+window_data_pts-1)); % Getting median power within window.
    end
    
    base_windows=(basetime-window)/winstep+1; % Total # windows in basetime.
    infuse_windows=(infusetime-window)/winstep+1;
    total_windows=length(InfuseBeta);
    
    Basebeta=nanmedian(InfuseBeta(1:base_windows)); % This is scalar?
    
    %% Finding changes in beta power, looking in sliding windows of length timetottest.
    timetottest=50; % 50sec to compare to baseline
    
    Betachangettest=nan(1,totalpoints-timetottest/winstep);
    Meanbeta=nan(1,totalpoints-timetottest/winstep);
    incbeta=[];
    decbeta=[];
    for i=1:(totalpoints-timetottest/winstep)
        Betachangettest(i)=ttest(InfuseBeta(i:(i+timetottest/winstep-1)), Basebeta, 0.05, 'both'); % p=0.05 or 0.01; most cases is 0.01
        Meanbeta(i)=nanmean(InfuseBeta(i:(i+timetottest/winstep-1)));
        % If mean beta is bigger than baseline beta, and the two-sided t-test comes back significant, there's a beta increase.
        if Meanbeta(i)>Basebeta && Betachangettest(i)==1
            incbeta(end+1)=i;
        end
        % If mean beta is smaller than baseline beta, and the two-sided t-test comes back significant, there's a beta increase.
        if Meanbeta(i)<Basebeta && Betachangettest(i)==1
            decbeta(end+1)=i;
        end
    end
    
    %% To eliminate false identified changes.
    for i=1:2
        
        if i==1
            tempbeta=incbeta;
        else
            tempbeta=decbeta;
        end
        
        % Creating vector of time indices where beta is changed.
        Tempchange=zeros(1,length(Betachangettest));
        Tempchange(tempbeta)=1;
        
        stringent=[1 1 1 1 1 1 1 1]; %%5 consecutive windows, that is 8*winstep=40sec.
        
        % Convolving index of change indices with square kernel (i.e., taking time avg. over 8 windows).
        TempBeta=conv(Tempchange, stringent, 'valid');
        TempBeta=horzcat(zeros(1, length(stringent)-1), TempBeta);
        
        TrueChange{i}=find(TempBeta>length(stringent)-2); % Want change to be over more than 6*winstep = 30 sec.
        if isempty (TrueChange{i})
            BetaOutPut(i,1)=NaN;  %latency of beta change onset
            BetaOutPut(i,2)=NaN;  %peak of beta power change
        else
            BetaOutPut(i,1)=TrueChange{i}(1)*winstep; % This is a scalar
            if i==1
                BetaOutPut(i,2)=max(Meanbeta(TrueChange{i}));
            end
            if i==2
                BetaOutPut(i,2)=min(Meanbeta(TrueChange{i}));
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
