function PD_decimate(prefix,plot_opt)

load([prefix, '_all_data.mat'])

PD_dec = decimate(PD_data, 10);

PD_dec = decimate(PD_dec, 2);

PD_dec = detrend(PD_dec, 'linear');

sampling_freq = sampling_rate/20;

hanning_window = hanning(20*sampling_freq);
hanning_window = hanning_window/sum(hanning_window);

PD_dec_reflected = [flipud(PD_dec(1:5*sampling_rate)); PD_dec; fliplr(PD_dec(1:5*sampling_rate))];

PD_trend = conv(PD_dec_reflected,hanning_window,'same');

PD_trend = PD_trend((5*sampling_rate + 1):(end - 5*sampling_rate));

if plot_opt > 0 
    
    t = (1:length(PD_dec))/(60*sampling_freq);
    
    figure;
    
    subplot(2,1,1)
    
    plot(t,PD_dec,'k',t,PD_trend,'r')
    
    legend({'Data','Trend'})
    
    title([prefix,', Data & Trend'])
    
    xlabel('Time (min.)')
    
end

PD_dec = PD_dec - PD_trend;

save([prefix,'_all_data_dec.mat'],'PD_dec','sampling_freq')

if plot_opt > 0
   
    subplot(2,1,2)
    
    plot(t,PD_dec)
    
    title([prefix,', All Data (Decimated & Detrended)'])
    
    xlabel('Time (min.)')
    
    save_as_pdf(gcf,[prefix,'_all_data_dec'])
    
end