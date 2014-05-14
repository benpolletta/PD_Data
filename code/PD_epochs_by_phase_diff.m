function PD_epochs_by_phase_diff(phase_lo,phase_hi,segment_length)

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    phase_name = [folder,'/',prefix,'_all_channel_data_dec_phase',num2str(phase_lo),'to',num2str(phase_hi)];
    
    load([folder,'/',prefix,'_all_channel_data_dec.mat'])
    
    signal_length = size(PD_dec,1);
    
    load([folder,'/',prefix,'_all_channel_data_dec_P_diff.mat'],'P_diff')
    
    phase_indicator = P_diff(:,3) > phase_lo & P_diff(:,3) < phase_hi;
    
    pi_smooth = conv(single(phase_indicator),ones(segment_length*sampling_freq,1)/(segment_length*sampling_freq),'same');
    
    [~,pi_peak_indicator] = spaced_peaks(pi_smooth,segment_length*sampling_freq,0);
    
    phase_indices = pi_peak_indicator & (pi_smooth >= .5);% & (pi_smooth >= quantile(pi_smooth,.5));
    
    save([phase_name,'.mat'],'phase_indicator','pi_smooth','phase_indices')
    
    epoch_list_fid = fopen([phase_name,'_epochs.list'],'w');
    
    epoch_phase_fid = fopen([phase_name,'_pi_smooth.txt'],'w');
    
    wincenters = find(phase_indices == 1);
    
    no_windows = length(wincenters);
    
    for w = 1:no_windows
        
        win_center = wincenters(w);
        
        fprintf(epoch_phase_fid,'%f\n',pi_smooth(win_center));
        
        segment_start = max(1,win_center-floor(segment_length*sampling_freq/2));
    
        segment_end = min(signal_length,win_center+floor(segment_length*sampling_freq/2));
    
        epoch_data = PD_dec(segment_start:segment_end,:);

        epoch_name = [phase_name,'_epoch',num2str(w),'.txt'];
        
        fid = fopen(epoch_name,'w');
        
        fprintf(fid,'%f\t%f\n',epoch_data');
        
        fclose(fid);
        
        fprintf(epoch_list_fid,'%s\n',epoch_name);
        
    end

    fclose(epoch_phase_fid);
    
    fclose(epoch_list_fid);
    
end