function PD_plot_wav(subjects_mat, epoch_length)

load(subjects_mat)

sampling_freq = 1000;

freqs = 1:200;

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

[r, c] = subplot_size(no_bands);

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

norm_labels = {'', ', % Total Power', ', % Baseline Power', ', Baseline Normalized % Total Power'};

stat_labels = {'Median', 'Mean'}; no_stats = length(stat_labels);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
   
    load([subj_name, '_', num2str(epoch_length/sampling_freq), 's_dec_wav.mat'])
    
    for n = 1:no_norms
        
        for s = 1:no_stats
            
            for ch = 1:2
            
                %% Plotting Spectrograms.
                
                figure((n - 1)*no_stats + s)
                
                subplot(2, 1, ch)
                
                imagesc(t_dec, freqs, Spec_dec(:, :, ch, n, s)')
                
                cl = caxis; caxis([cl(1) .25*cl(2)])
                
                axis xy
                
                xlabel('Time (s)'); ylabel('Hz');
                
                title(['Gabor Spectrogram of ', folder, ', ', chan_labels{ch}, ', ', stat_labels{s}, norm_labels{n}])
                
                grid on
                
            end
        
            try save_as_pdf(gcf, [subj_name, '_', num2str(epoch_length/sampling_freq), 's_wav', norms{n}, '_', stat_labels{s}]), end
            
            %% Plotting Band Power.
            
            figure(no_norms*no_stats + (n - 1)*no_stats + s)
            
            for b = 1:no_bands
                
                subplot(r, c, b)
                
                plot(t_dec, reshape(BP_dec(:, b, :, n, s), size(BP_dec, 1), size(BP_dec, 3)))
                
                axis tight
                
                xlabel('Time (s)')
                
                ylabel([stat_labels{s}, ' Power', norm_labels{n}])
                
                title(['Power, ', num2str(bands(b, 1)), ' - ', num2str(bands(b, 2)), ' Hz'])
                
                legend(chan_labels)
                
            end
            
            try save_as_pdf(gcf, [subj_name, '_', num2str(epoch_length/sampling_freq), 's_wav_BP', norms{n}, '_', stat_labels{s}]), end
            
        end
        
    end
    
end