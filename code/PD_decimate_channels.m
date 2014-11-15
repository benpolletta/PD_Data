function PD_decimate_channels(prefix,plot_opt)

load([prefix,'_all_channel_data.mat'])

figure;

for ch = 1:2
    
    PD_dec_temp = decimate(PD_data(:,ch), 10);
    
    PD_dec_temp = decimate(PD_dec_temp, 2);
    
    PD_dec_temp = detrend(PD_dec_temp, 'linear');
    
    sampling_freq = sampling_rate/20;
    
    hanning_window = hanning(20*sampling_freq);
    hanning_window = hanning_window/sum(hanning_window);
    
    PD_dec_temp_reflected = [flipud(PD_dec_temp(1:5*sampling_rate)); PD_dec_temp; fliplr(PD_dec_temp(1:5*sampling_rate))];
    
    PD_trend = conv(PD_dec_temp_reflected,hanning_window,'same');
    
    PD_trend = PD_trend((5*sampling_rate + 1):(end - 5*sampling_rate));
    
    if plot_opt > 0
        
        t = (1:length(PD_dec_temp))/(60*sampling_freq);
        
        subplot(3,1,ch)
        
        plot(t,PD_dec_temp,'k',t,PD_trend,'r')
        
        legend({'Data','Trend'})
        
        title([prefix,', Data & Trend'])
        
        xlabel('Time (min.)')
        
    end
    
    PD_dec_temp = PD_dec_temp - PD_trend;
    
    PD_dec(:,ch) = PD_dec_temp;
    
end

save([prefix,'_all_channel_data_dec.mat'],'PD_dec','sampling_freq')

if plot_opt > 0
   
    subplot(3,1,3)
    
    plot(t,PD_dec)
    
    title([prefix,', All Data (Decimated & Detrended)'])
    
    legend({'Channel 1','Channel 2'})
    
    xlabel('Time (min.)')
    
    save_as_pdf(gcf,[prefix,'_all_channel_data_dec'])
    
end