function make_6OHDA_spectrogram_figures(subject_no, subject_name, subject_prefix)
    
group_prefix = 'st_m1_6OHDA';
data = load([subject_name, '/', subject_prefix, '_wt.mat']);
Spec = data.Spec;

for ch = 1:2
    
    figure
    
    %% Plotting spectrogram for whole data recording.
    
    Ch_spec = abs(Spec(:, :, ch));
    
    dims = size(Ch_spec);
    
    dec_factor = 1; s_rate = 500/dec_factor;
    
    Spec_dec = nanmean(reshape(Ch_spec', dims(2), dec_factor, dims(1)/dec_factor), 2);
    Spec_dec = permute(Spec_dec, [1 3 2]);
    
    t = ((1:length(Spec_dec))/s_rate)/60;
    
    Spec_dec = Spec_dec(1:80, :);
    
    subplot(4, 1, 1)
    
    imagesc(t, 1:80, Spec_dec(1:80, :))
    
    title([subject_name, ', Channel ', num2str(ch)])
    
    ylabel('Frequency (Hz)')
    
    caxis([0 10])
    
    axis xy
    
    hold on
    
    plot([0 0], [0 80], ':w', 'LineWidth', 2)
    
    plot([t(1) t(end)], [30 30], ':w', 'LineWidth', 2)
    
    plot([t(1) t(end)], [8 8], ':w', 'LineWidth', 2)
    
    load([group_prefix, '_pct_BP_high_2.5_min_secs.mat'])
    
    starts = All_bp_max_start(subject_no, :, 3);
    start = starts(ch);
    start = (start/s_rate)/60;
    
    ends = All_bp_max_end(subject_no, :, 3);
    ending = ends(ch);
    ending = (ending/s_rate)/60;
    
    plot([1 1]*start, [0 80], 'y')
    plot([1 1]*ending, [0 80], 'y')
    
    %% Plotting LFP for whole data recording.
    
    load([subject_name, '/', subject_prefix, '_all_channel_data_dec.mat']);
    
    subplot(4, 1, 2)
    
    plot(t, PD_dec(:, ch), 'k')
    
    box off
    
    axis tight
    
    xlabel('Time (min.) Rel. Infusion')
    
    ylabel('LFP (mV)')
    
    %% Plotting spectrogram & LFP for 10 seconds.
    
    midpoint = (start + ending)/2;
    
    start_time = midpoint - 5/60; end_time = midpoint + 5/60;
    
    t_index = t >= start_time & t <= end_time;
    
    t_interval = t(t_index)*60;
    
    subplot(4, 1, 3)
    
    imagesc(t_interval, 1:80, Ch_spec(t_index, 1:80)')
    
    title('Densest Beta')
    
    ylabel('Frequency (Hz)')
    
    caxis([0 10])
    
    axis xy
    
    hold on
    
    plot([0 0], [0 80], ':w', 'LineWidth', 2)
    
    plot([t_interval(1) t_interval(end)]*60, [30 30], ':w', 'LineWidth', 2)
    
    plot([t_interval(1) t_interval(end)]*60, [8 8], ':w', 'LineWidth', 2)
    
    subplot(4, 1, 4)
    
    plot(t_interval, PD_dec(t_index, ch), 'k')
    
    box off
    
    axis tight
    
    xlabel('Time (sec.) Rel. Infusion')
    
    ylabel('LFP (mV)')
    
    try
        
        save_as_eps(gcf, [subject_name, '/', subject_name, '_spec_for_paper_ch', num2str(ch)])
        
    catch
        
        save(gcf, [subject_name, '/', subject_name, '_spec_for_paper_ch', num2str(ch), '.fig'])
        
    end
    
end
    
end
