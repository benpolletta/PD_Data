function PD_plot_pmtm(subjects_mat, epoch_secs, window_secs)

close('all')

load(subjects_mat)

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

[r, c] = subplot_size(no_bands);

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

norm_labels = {'', ', % Total Power', ', % Baseline Power', ', Baseline Normalized % Total Power'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
   
    All_data = load([subj_name, '_', num2str(epoch_secs), 's_epoch_pmtm.mat']);
    
    load([subj_name, '_', num2str(epoch_secs), 's_pmtm_', num2str(window_secs), 's_stats.mat']);
    
    lower_test = p_less <= .01/(2*sum_all_dimensions(~isnan(p_less))); lower_test = +lower_test;
    
    lower_test(lower_test == 0) = nan;
    
    upper_test = p_greater <= .05/(2*sum_all_dimensions(~isnan(p_greater))); upper_test = +upper_test;
    
    upper_test(upper_test == 0) = nan;
    
    t = All_data.t;
    
    freqs = All_data.freqs;
    
    plot_freq_indices = freqs >= 0 & freqs <= 200;
    
    for n = 1:no_norms
            
        %% Plotting Spectrograms.
        
        Spec_data = getfield(All_data, ['Spec', norms{n}]);
        
        for ch = 1:2
            
            figure((fo - 1)*no_norms*2 + n)
            
            subplot(2, 1, ch)
            
            if n < 3
                
                Spec_plot = zscore(Spec_data(:, plot_freq_indices, ch))';
                
            else
                
                Spec_plot = Spec_data(:, plot_freq_indices, ch)';
                
            end
            
            imagesc(t, freqs(plot_freq_indices), Spec_plot)
            
            cl = caxis; caxis([cl(1) .25*cl(2)])
            
            axis xy
            
            xlabel('Time (s)'); ylabel('Hz');
            
            title(['Multi-taper Spectrogram of ', folder, ', ', chan_labels{ch}, norm_labels{n}])
            
            grid on
            
        end
        
        try, save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_epochs_pmtm', norms{n}]), end
        
        %% Plotting Band Power.
        
        BP_data = getfield(All_data, ['BP', norms{n}]);
        
        figure((fo - 1)*no_norms*2 + no_norms + n)
        
        for b = 1:no_bands
            
            handle = subplot(r, c, b);
            
            plot(t, zscore(reshape(BP_data(:, b, :), size(BP_data, 1), 2)))
            
            add_stars(handle, t_win, [lower_test(:, :, b, n) upper_test(:, :, b, n)], [0 0 1 1], [0 0 1; 0 .5 0; 0 0 1; 0 .5 0])
            
            xlabel('Time (s)')
            
            ylabel(['Power', norm_labels{n}])
            
            title([folder, ', ', num2str(bands(b, 1)), ' - ', num2str(bands(b, 2)), ' Hz'])
            
            if b == 1
            
                legend(chan_labels, 'Location', 'NorthWest')
            
            end
            
        end
        
        try, save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_epochs_' , num2str(window_secs), 's_stats_pmtm_BP', norms{n}]), end
        
    end
    
    % close('all')
    
end