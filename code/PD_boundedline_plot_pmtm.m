function PD_boundedline_plot_pmtm(subjects_mat, epoch_secs, window_secs)

close('all')

load(subjects_mat), no_subjects = length(folders);

load([folders{1}, '/', prefixes{1}, '_', num2str(epoch_secs), 's_epoch_pmtm.mat'], 'freqs')

plot_freq_indices = freqs >= 0 & freqs <= 200; plot_freqs = freqs(plot_freq_indices);

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

norm_labels = {'', ', % Total Power', ', % Baseline Power', ', Baseline Normalized % Total Power'};

epochs_per_window = window_secs/epoch_secs;

load([subjects_mat(1:(end - length('subjects.mat'))), '_boundedline_collect_pmtm.mat'])

for b = 3:3 % 1:no_bands
    
    band_indices = plot_freqs >= bands(b, 1) & plot_freqs <= bands(b, 2);
    
    for n = 1:no_norms
        
        for ch = 1:2
            
            figure
            
            subplot(2, 1, 1)
            
            mean_for_plot = nanmean(All_pow(:, :, [0 3] + ch, ch, b, n));
            
            mean_for_plot = reshape(mean_for_plot, size(mean_for_plot, 2), size(mean_for_plot, 3));
            
            se_for_plot = nanstd(All_pow(:, :, [0 3] + ch, ch, b, n))/sqrt(size(All_pow, 1));
            
            se_for_plot = reshape(se_for_plot, size(se_for_plot, 2), size(se_for_plot, 3));
            
            se_for_plot = prep_for_boundedline(se_for_plot(plot_freq_indices, :));
            
            boundedline(plot_freqs, mean_for_plot(plot_freq_indices, :), se_for_plot)
            
            legend({'Pre','Post'})
            
            axis tight
            
            title([sprintf('%.0f - %.0f Hz', bands(b, :)), norm_labels{n}, ', Minute of Highest ',...
                chan_labels{ch}, ' Power, St. Dev. Over ', num2str(epoch_secs), ' Epochs'])
            
            subplot(2, 1, 2)
            
            mean_for_plot = nanmean(All_mean(:, :, [0 3] + ch, ch, b, n));
            
            mean_for_plot = reshape(mean_for_plot, size(mean_for_plot, 2), size(mean_for_plot, 3));
            
            se_for_plot = nanstd(All_mean(:, :, [0 3] + ch, ch, b, n))/sqrt(size(All_mean, 1));
            
            se_for_plot = reshape(se_for_plot, size(se_for_plot, 2), size(se_for_plot, 3));
            
            se_for_plot = prep_for_boundedline(se_for_plot(plot_freq_indices, :));
            
            boundedline(plot_freqs, mean_for_plot(plot_freq_indices, :), se_for_plot)
            
            legend({'Pre','Post'})
            
            axis tight
            
            title([sprintf('%.0f - %.0f Hz', bands(b, :)), norm_labels{n}, ', Minute of Highest ',...
                chan_labels{ch}, ' Power, St. Dev. Over Subjects'])
            
%             for fo = 1:no_subjects
%                 
%                 folder = folders{fo};
%                 
%                 subplot(no_subjects + 1, 2, 2 + 2*fo - 1)
%                 
%                 All_pow((fo - 1)*epochs_per_window + (1:epoch_per_window), :, :, ch, b) = Spec_window;
%                 
%                 All_mean(fo, :, :, ch, b) = mean(Spec_window);
%                 
%             end
            
        end
        
    end
    
end
                
%                 figure
%             
%                 subplot(no_subjects, 2, 2*fo - 1)
%             
%             imagesc(t, freqs(plot_freq_indices), Spec_plot)
%             
%             cl = caxis; caxis([cl(1) median(max(Spec_plot, [], 2))])
%             
%             axis xy
%             
%             xlabel('Time (s)'); ylabel('Hz');
%             
%             title(['Multi-taper Spectrogram of ', folder, ', ', chan_labels{ch}, norm_labels{n}])
%             
%             grid on
%             
%         end
%         
%         try, save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_epochs_pmtm', norms{n}]), end
%         
%         %% Plotting Band Power.
%         
%         BP_data = getfield(All_data, ['BP', norms{n}]);
%         
%         BP_data(artifact_indicator, :, :) = nan;
%         
%         figure((fo - 1)*no_norms*2 + no_norms + n)
%             
%             handle = subplot(r, c, b);
%             
%             plot(t, nanzscore(reshape(BP_data(:, b, :), size(BP_data, 1), 2)))
%             
%             add_stars(handle, t_win, [lower_test(:, :, b, n) upper_test(:, :, b, n)], [0 0 1 1], [0 0 1; 0 .5 0; 0 0 1; 0 .5 0])
%             
%             xlabel('Time (s)')
%             
%             ylabel(['Power', norm_labels{n}])
%             
%             title([folder, ', ', num2str(bands(b, 1)), ' - ', num2str(bands(b, 2)), ' Hz'])
%             
%             if b == 1
%             
%                 legend(chan_labels, 'Location', 'NorthWest')
%             
%             end
%             
%         end
%         
%         try, save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_epochs_' , num2str(window_secs), 's_stats_pmtm_BP', norms{n}]), end
%         
%     end
%     
%     % close('all')
%     
% end
% 
% end
% 
% function times = get_times_pow_p(t_win_pow, power, t_win_p, p_less, p_greater)          
% 
%             [~, max_pre_index] = max(power(t_win_pow < 0, :, b, n);
%             
%             increase_pre = t_win_pow(max_pre_index);
%             
%             [~, min_pre_index] = min(power(t_win_pow < 0, :, b, n);
%             
%             decrease_pre = t_win_pow(min_pre_index);
%             
%             [~, increase_post_index] = max(p_greater(:, :, b, n));
%             
%             increase_post = t_win_p(increase_post_index);
%             
%             [~, decrease_post_index] = max(p_less(:, :, b, n));
%             
%             decrease_post = t_win_p(decrease_post_index);
%             
%             times = [increase_pre'; decrease_pre'; increase_post'; decrease_post'];
%             
% end