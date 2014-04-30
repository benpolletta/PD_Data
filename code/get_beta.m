present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

basetimes = [300 1200 1800 600 1800];

infusetimes = [390 240 300 450 510];

load('BetaTimes')

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    load([prefix,'_all_channel_data_dec.mat']);
   
    beta_start_index = BetaTimes(fo,1)*sampling_freq
    
    beta_end_index = min(length(PD_dec), BetaTimes(fo,2)*sampling_freq)
    
    beta_data = PD_dec(beta_start_index:beta_end_index,:);
    
    save([prefix,'_beta_data.mat'],'beta_data')
    
    [~,~,A,~] = filter_wavelet_Jan(beta_data(:,1),'bands',1:40,'sampling_freq',sampling_freq,'filename',prefix);
    
    figure;
    
    imagesc((1:length(beta_data))/sampling_freq,1:40,zscore(A)')
    
    axis xy
    
    title(prefix)
    
    xlabel('Time (s)')
    
    ylabel('Frequency (Hz)')
    
    saveas(gcf,[prefix,'_wav_1-40_Hz.fig'])
    
    beta_sig = eegfilt(beta_data',sampling_freq,15,30);
    
    H = hilbert(beta_sig);
    
    figure;
    
    plot((1:length(beta_sig))/sampling_freq,abs(H))
    
    legend({'Channel 1','Channel 2'})
    
    saveas(gcf,[prefix,'_bp_15-30_Hz.fig'])
    
    cd (present_dir)
    
end