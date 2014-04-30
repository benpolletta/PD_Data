function BetaOut = Park_pop_betapower_v3 (chData, basetime, infusetime)

    srate = 1000; % Downsampled sampling rate.
    binrate = 20; % How many timebins per second.
    timebinlength = srate/binrate; % If 10, 100ms per point; if 20, 50ms each points; if 50, 20ms each point.
    betaPower = (chData.beta).^2; % Squared magnitude of FFT; should be sampled at srate.
    totaltime = length(chData.beta)/binrate; % Must be one measurement of beta power per timebin?
    
    %%
    window = 5; % In sec.
    winstep = 5; % In sec.
    
    winstep_data_pts = winstep*binrate;
    window_data_pts = window*binrate;
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
    
    Betachangettest=nan(1,total_windows-timetottest/winstep);
    Meanbeta=nan(1,total_windows-timetottest/winstep);
    incbeta=[];
    decbeta=[];
    for i=1:(total_windows-timetottest/winstep)
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
        
        TrueChange{i}=TempBeta>length(stringent)-2; % Want change to be over more than 6*winstep = 30 sec.
        if ~any(TrueChange{i})>0
            BetaOutPut(i,1)=NaN;  % Latency of beta change onset.
            BetaOutPut(i,2)=NaN;  % Peak of beta power change.
        else
            
            Change_diff = diff(TrueChange{i});
            Change_on = find(Change_diff == 1);
            Change_off = find(Change_diff == -1);
            
            Temp_meanbeta = Meanbeta;
            Temp_meanbeta(TrueChange{i} == 0) = 0;
            
            if i==1
                [~, BetaMaxIndex] = max(Temp_meanbeta);
            else
                [~, BetaMaxIndex] = min(Temp_meanbeta);
            end
                
            Change_on_index = max(Change_on(Change_on <= BetaMaxIndex));
            Change_off_index = min(Change_off(Change_off >= BetaMaxIndex));
            
            BetaOutPut(i,1) = Change_on_index*winstep-8*winstep; % This is a scalar.
            BetaOutPut(i,2) = Change_off_index*winstep-winstep;
            % This is max & min of beta power, but I don't need this.
            % if i==1
            %     BetaOutPut(i,2)=max(Meanbeta(TrueChange{i})); % But this seems to be the max. beta power.
            % end
            % if i==2
            %     BetaOutPut(i,2)=min(Meanbeta(TrueChange{i})); % And this seems to be the min. beta power.
            % end
        
        end
        
    end
    
    recoverbeta=nanmedian(InfuseBeta(total_windows-base_windows+1:total_windows));
    
    BetaOut = [BetaOutPut(1,:) BetaOutPut(2,:)];
    
    % BetaOut=[Basebeta BetaOutPut(1,:) BetaOutPut(2,:) recoverbeta]; 
    % So, normally, this gives mean beta power during baseline, max mean
    % beta power (over a 50 s window) during beta increase, min mean beta
    % power (over a 50 s window) during beta increase, max etc. etc. during
    % beta decrease, min etc. etc. during beta decrease, and finally mean
    % beta power over all not baseline times. I don't know why you would
    % need this last thing. I just want the latencies.
    
    % I don't know what the fuck this is.
%     BetaControl=[Basebeta 750 nanmedian(InfuseBeta(150:150+base_windows-1)) NaN NaN recoverbeta];
    
    figure;
    plot(1:winstep:length(InfuseBeta)*winstep, InfuseBeta, 'r-')
    hold on;
    scatter(1:winstep:length(Betachangettest)*winstep, Betachangettest*(max(InfuseBeta)-min(InfuseBeta))+min(InfuseBeta))
    plot(repmat(BetaOutPut(1,:),2,1), repmat([min(InfuseBeta); max(InfuseBeta)],1,2),'k')
    plot(repmat(BetaOutPut(2,:),2,1), repmat([min(InfuseBeta); max(InfuseBeta)],1,2),'g')
    
%     figure (1)
%     x=0:timebinlength/1000:length(chData.beta)*timebinlength/1000;
%     y=1:200;
%     imagesc(x, y, log(chData.allfreq))
%     axis xy
%     axis ([0 length(chData.beta)*timebinlength/1000 4 70])
%     xlabel ('time (s)');
%     ylabel ('Frequency (Hz)')
%     colorbar;
     
end
