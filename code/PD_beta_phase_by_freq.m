function PD_beta_phase_by_freq(no_p_bins,f_bins)

sampling_freq = 1000;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

channel_label = {'Striatum','Motor Ctx.'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    phase_name = [folder,'/',prefix,'_beta_data_PbyF'];
    
    load([folder,'/',prefix,'_beta_data.mat'])
    
    BP = eegfilt(beta_data',sampling_freq,10,30);
    
    H = hilbert(BP'); A = abs(H); P = angle(H);
    
    P_unwrapped = unwrap(P);
    
    P_diff = angle(exp(sqrt(-1)*(-diff(P_unwrapped,[],2))));
    
    F = diff(P_unwrapped)/(2*pi*(1/sampling_freq));
    
    F_smooth = nan(size(F)); P_diff_smooth = nan(size(P));
    
    for ch = 1:2
       
        F_smooth(:,ch) = conv(F(:,ch),ones(50,1)/50,'same');
        
    end
       
    P_diff_smooth = angle(conv(exp(sqrt(-1)*P_diff),ones(500,1)/500,'same'))/pi;
    
    figure;
    
    subplot(2,1,1)
    
    plot((1:size(F,1))'/sampling_freq,F_smooth)
    
    axis tight
    
    title(folder)
    
    xlabel('Time (s)')
    
    ylabel('Frequency (Hz)')
    
    legend({'Striatal','Motor'})
    
    subplot(2,1,2)
    
    plot((1:size(P,1))'/sampling_freq,P_diff_smooth)
    
    axis tight
    
    xlabel('Time (s)')
    
    ylabel('Phase Difference (Striatal - Motor, \pi)')
    
    MR_vec = nan(2, length(f_bins)-1);
    
    p_hist = nan(no_p_bins+1, length(f_bins)-1, 2);
    
    freqs = nan(2, length(f_bins)-1);
    
    figure;
    
    for ch = 1:2
        
        subplot(1,2,ch)
        
        [MR_vec(ch,:), p_hist(:,:,ch), freqs(ch,:)] = rose_plot(P_diff, F_smooth(:,ch), no_p_bins, f_bins);
        
        title({[folder,', Phase Difference (Striatal - Motor)'];['By ',channel_label{ch},' Frequency']})
        
    end
    
    save([phase_name,'.mat'],'H','A','P','P_diff','F','F_smooth','MR_vec','p_hist','freqs')
    
    save_as_pdf(gcf,phase_name)
    
end