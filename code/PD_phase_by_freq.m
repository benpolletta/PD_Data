function PD_phase_by_freq(no_p_bins,f_bins)

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    phase_name = [folder,'/',prefix,'_all_channel_data_dec_PbyF'];
    
    load([folder,'/',prefix,'_all_channel_data_dec_HAP.mat'],'P')
    
    F = diff(unwrap(P(:,:,3)))/(2*pi*(1/1000));
    
    F_smooth = nan(size(F));
    
    for ch = 1:2
       
        F_smooth(:,ch) = conv(F(:,ch),ones(50,1)/50,'same');
        
    end
    
    figure;
    
    plot((1:size(F,1))'/1000,F_smooth)
    
    load([folder,'/',prefix,'_all_channel_data_dec_P_diff.mat'],'P_diff')
    
    MR_vec = nan(2, length(f_bins)-1);
    
    p_hist = nan(no_p_bins+1, length(f_bins)-1, 2);
    
    freqs = nan(2, length(f_bins)-1);
    
    figure;
    
    for ch = 1:2
        
        subplot(1,2,ch)
        
        [MR_vec(ch,:), p_hist(:,:,ch), freqs(ch,:)] = rose_plot(P_diff, F_smooth(:,ch), no_p_bins, f_bins);
        
    end
    
    save([phase_name,'.mat'],'F','F_smooth','MR_vec','p_hist','freqs')
    
    save_as_pdf(gcf,phase_name)
    
end