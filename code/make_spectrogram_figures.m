function make_spectrogram_figures(subject_no, subject_name, basetime, striatum_channel)
    
data = load([subject_name, '/', subject_name([1 2 4:end]), '_wt.mat']);
Spec = data.Spec;

for i = 1:2
    
    figure
    
    %% Plotting spectrogram for whole data recording.
    
    Str_spec = abs(Spec(:, :, i));
    
    dims = size(Str_spec);
    
    dec_factor = 1; s_rate = 500/dec_factor;
    
    Spec_dec = nanmean(reshape(Str_spec', dims(2), dec_factor, dims(1)/dec_factor), 2);
    Spec_dec = permute(Spec_dec, [1 3 2]);
    
    t = ((1:length(Spec_dec))/s_rate - basetime)/60;
    
    Spec_dec = Spec_dec(1:80, t <= 30);
    t(t > 30) = [];
    
    subplot(4, 1, 1)
    
    imagesc(t, 1:80, Spec_dec(1:80, :))
    
    title([subject_name, ', Channel ', num2str(i)])
    
    ylabel('Frequency (Hz)')
    
    caxis([0 10])
    
    axis xy
    
    hold on
    
    plot([0 0], [0 80], ':w', 'LineWidth', 2)
    
    plot([t(1) t(end)], [30 30], ':w', 'LineWidth', 2)
    
    plot([t(1) t(end)], [8 8], ':w', 'LineWidth', 2)
    
    load('st_m1_pct_BP_high_2.5_min_secs.mat')
    
    starts = reshape(All_bp_max_start(subject_no, :, 3, :), 2, 2);
    starts = starts(striatum_channel, :);
    
    ends = reshape(All_bp_max_end(subject_no, :, 3, :), 2, 2);
    ends = ends(striatum_channel, :);
    
    colors = {'g', 'r'};
    
    for i = 1:2
        
        plot(([1 1]*starts(i)/500 - basetime)/60, [0 80], colors{i})
        plot(([1 1]*ends(i)/500 - basetime)/60, [0 80], colors{i})
        
    end
    
    %% Plotting LFP for whole data recording.
    
    load([subject_name, '/', subject_name([1 2 4:end]), '_all_channel_data_dec.mat']);
    
    subplot(4, 1, 2)
    
    plot(t, PD_dec(t <= 30, i), 'k')
    
    box off
    
    axis tight
    
    xlabel('Time (min.) Rel. Infusion')
    
    ylabel('LFP (mV)')
    
    %% Plotting spectrogram & LFP for 10 seconds.
    
    pd_labels = {'Pre-Infusion', 'Post-Infusion'};
    
    for i = 1:2
        
        midpoint = ((starts(i)/500 - basetime)/60 + (ends(i)/500 - basetime)/60)/2;
        
        start_time = midpoint - 5/60; end_time = midpoint + 5/60;
        
        t_index = t >= start_time & t <= end_time;
        
        t_interval = t(t_index);
        
        subplot(4, 2, 4 + i)
        
        imagesc(t_interval, 1:80, Str_spec(t_index, 1:80)')
        
        title(pd_labels{i})
        
        ylabel('Frequency (Hz)')
        
        caxis([0 10])
        
        axis xy
        
        hold on
        
        plot([0 0], [0 80], ':w', 'LineWidth', 2)
        
        plot([t_interval(1) t_interval(end)], [30 30], ':w', 'LineWidth', 2)
        
        plot([t_interval(1) t_interval(end)], [8 8], ':w', 'LineWidth', 2)
        
        subplot(4, 2, 6 + i)
        
        plot(t_interval, PD_dec(t_index, i), 'k')
        
        box off
        
        axis tight
        
        xlabel('Time (min.) Rel. Infusion')
        
        ylabel('LFP (mV)')
        
        try
            
            save_as_pdf(gcf, [subject_name, '/', subject_name, '_spec_for_paper_ch', num2str(i)])
            
        catch
            
            saveas(gcf, [subject_name, '/', subject_name, '_spec_for_paper_ch', num2str(i), '.fig'])
            
        end
        
    end
    
end
