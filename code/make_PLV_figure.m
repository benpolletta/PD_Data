function make_PLV_figure

load(['STR_M1_subjects.mat'], 'chan_labels', 'pd_labels', 'folders')

group_flags = {'M1_not_increased', 'M1_increased'};

group_titles = {'M1 \beta >'; 'M1 \beta \leq'};

channel_prefixes = {'STR_w_M1', 'M1'};

no_chans = length(chan_labels);

peak_suffix = '_kmeans';

no_pds_plotted = 2;

figure

freq_limit = 200;

p_val = .01;

for group = 1:2

    %% Plotting magnitude of coherence.

    load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_', group_flags{group} '_Coh_sec_pct_data_for_plot.mat'])
    
    PLV_mean = All_mean_mean; PLV_ci = All_mean_ci;
    
    subplot(3, 3, 3*(group - 1) + 3) % 6 + group)
    
    boundedline((1:freq_limit)', PLV_mean(1:freq_limit, 1:no_pds_plotted), prep_for_boundedline(PLV_ci(1:freq_limit, 1:no_pds_plotted)))
    
    axis tight
    
    y_lims(group, :) = ylim;
    
end

y_extremes(3, :) = [min(y_lims(:, 1)) max(y_lims(:, 2))];

for group = 1:2
    
    load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_', group_flags{group} '_Coh_sec_pct_data_for_plot.mat'])
    
    PLV_mean = All_mean_mean; PLV_ci = All_mean_ci;
    
    [sig_lower, sig_higher] = find_sig(PLV_mean(:, 1:2), PLV_ci(:, 1:2));
    
    subplot(3, 3, 3*(group - 1) + 3) % 6 + group)
    
    ylim(y_extremes(3, :))
    
    add_stars(gca, (1:freq_limit)', logical(sig_lower(1:freq_limit)), 0, [1 .5 0])
    
    add_stars(gca, (1:freq_limit)', logical(sig_higher(1:freq_limit)), 1, [1 0 0])
    
    hold on
    
    % plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')
    
    % plot((1:freq_limit)', zeros(1, freq_limit), ':k')
    
    % axis tight
    
    set(gca, 'FontSize', 16)
    
    xlabel('Freq. (Hz)', 'FontSize', 16)
    
    if group == 1
    
        title({'Phase-Locking Magnitude'; '(%\Delta BL, Mean \pm 95% CI)'}, 'FontSize', 16)
    
    end
    
end

%% Plotting spectra.

for ch = 1:2
    
    for group = 1:2
        
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
        
        All_mean_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
        
        subplot(3, 3, 3*(group - 1) + ch) % 3*(ch - 1) + group)
        
        boundedline((1:freq_limit)', All_mean_mean(1:freq_limit, :), prep_for_boundedline(All_mean_ci(1:freq_limit, :)))
        
        axis tight
        
        y_lims(group, :) = ylim;
        
    end

    y_extremes(ch, :) = [min(y_lims(:, 1)) max(y_lims(:, 2))];
    
    for group = 1:2
        
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
        
        All_mean_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
        
        [sig_lower, sig_higher] = find_sig(All_mean_mean, All_mean_ci);
        
        subplot(3, 3, 3*(group - 1) + ch) % 3*(ch - 1) + group)
        
        ylim(y_extremes(ch, :))
        
        add_stars(gca, (1:freq_limit)', logical(sig_lower(1:freq_limit)), 0, [1 .5 0])
        
        add_stars(gca, (1:freq_limit)', logical(sig_higher(1:freq_limit)), 1, [1 0 0])
        
        % plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')
        
        set(gca, 'FontSize', 16)
        
        % xlabel('Freq. (Hz)', 'FontSize', 16)
        
        if group == 1
            
            title({[chan_labels{ch}, ' Power']; '(% \Delta BL, Mean \pm 95% CI)'}, 'FontSize', 16)
            
        end
        
        if ch == 1
            
            y = ylabel(group_titles{group}, 'FontSize', 20, 'Rotation', 0); % title(chan_labels{1}, 'FontSize', 20)
            
            set(y, 'Units', 'Normalized', 'Position', [-0.35 0.4 0])
            
        end
        
    end
    
end

%% Computing stats between post-infusion increases.

load('M1_groups.mat')

for group = 1:2
    
    no_excluded = length(M1_groups{group});
    
    folder_index{group} = ones(1, length(folders));
    
    for e = 1:no_excluded
    
        folder_index{group} = folder_index{group} - strcmp(folders, M1_groups{group}{e});
        
    end

end

for ch = 1:2
    
    load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat'])
    
    p_vals = nan(freq_limit, 2);
    
    for f = 1:freq_limit
        
        [~, p_vals(f, 1)] = ttest2(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'right');
        
        [~, p_vals(f, 2)] = ttest2(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'left');
        
    end
    
    test = p_vals < p_val;
    
    mean_mat = nan(freq_limit, 2);
    
    for group = 1:2
        
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
    
        mean_mat(:, group) = All_mean_mean(:, 2);
        
    end
        
    subplot(3, 3, 3*2 + ch) % subplot(3, 3, 3*(ch - 1) + 3)
    
    plot((1:freq_limit)', -diff(mean_mat(1:freq_limit, :), [], 2), 'k', 'LineWidth', 1.5)
    
    ylim(y_extremes(ch, :))
    
    add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 1)), 0, [0 1 1])
    
    add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 2)), 1, [1 0 1])

    if ~any(test ~= 0), hold on, end
    
    plot((1:200)', zeros(200, 1), 'k:')
    
    if ch == 1
        
        y = ylabel({group_titles{1}; [' - ', group_titles{2}]}, 'FontSize', 20, 'Rotation', 0);
    
        set(y, 'Units', 'Normalized', 'Position', [-0.4 0.4 0])
        
    end
    
    set(gca, 'FontSize', 16)
    
    box off
    
end

load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_Coh_sec_pct_data_for_plot.mat'])

p_vals = nan(freq_limit, 2);

for f = 1:freq_limit
    
    [~, p_vals(f, 1)] = ttest2(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'right');
    
    [~, p_vals(f, 2)] = ttest2(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'left');
    
end

test = p_vals < p_val;

mean_mat = nan(freq_limit, 2);

for group = 1:2
    
    load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_', group_flags{group}, '_Coh_sec_pct_data_for_plot.mat'])
    
    mean_mat(:, group) = All_mean_mean(:, 2);
    
end

subplot(3, 3, 3*2 + 3)

plot((1:freq_limit)', -diff(mean_mat(1:freq_limit, :), [], 2), 'k', 'LineWidth', 1.5)

ylim(y_extremes(3, :))

add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 1)), 0, [0 1 1])

add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 2)), 1, [1 0 1])

if ~any(test ~= 0), hold on, end 

plot((1:200)', zeros(200, 1), 'k:')

set(gca, 'FontSize', 16)

box off

xlabel('Freq. (Hz)', 'FontSize', 16)

% %% Comparing striatum & motor cortex.
% 
% freq_limit = 200;
% 
% for b = 4
%     
%     All_mean_channels = nan(200, length(folder_index{1}), 2);
%     
%     for ch = 1:2
%         
%         load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands_kmeans_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat'])
%         
%         All_mean_channels(:, :, ch) = All_mean(:, :, 2);
%         
%     end
%     
%     for g = 1:2
%         
%         mean_mat = nan(freq_limit, 2);
%         
%         for ch = 1:2
%             
%             load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands_kmeans_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_', group_flags{g}, '_ch1_data_for_plot.mat'])
%             
%             mean_mat(:, ch) = All_mean_mean(:, 2);
%             
%         end
%         
%         subplot(3, 4, 4*(g - 1) + 4)
%         
%         plot((1:freq_limit)', -diff(mean_mat(1:freq_limit, :), [], 2), 'k', 'LineWidth', 1.5)
%         
%         ylim([min(y_extremes(1:2, 1)) max(y_extremes(1:2, 2))])
%         
%         % axis tight
%          
% %         y_lims(g, :) = ylim;
% %         
% %     end
% %     
% %     y_extremes = [min(y_lims(:, 1)) max(y_lims(:, 2))];
% %     
% %     for g = 1:2
%         
%         p_vals = nan(freq_limit, 2);
%         
%         for f = 1:freq_limit
%             
%             [~, p_vals(f, 1)] = ttest2(All_mean_channels(f, logical(folder_index{g}), 1)', All_mean(f, logical(folder_index{g}), 2)', 'tail', 'left');
%             
%             [~, p_vals(f, 2)] = ttest2(All_mean_channels(f, logical(folder_index{g}), 1)', All_mean(f, logical(folder_index{g}), 2)', 'tail', 'right');
%             
%         end
%         
%         test = p_vals < p_val;
%         
%         % h = subplot(3, 4, 4*(g - 1) + 4);
%         
%         % set(h, 'ylim', y_extremes)
%         
%         add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 1)), 0, [1 1 0])
%         
%         add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 2)), 1, [.5 0 .5])
%         
%         if ~any(test ~= 0), hold on, end
%         
%         plot((1:200)', zeros(200, 1), 'k:')
%         
%         if g == 1
%             
%             title({[' Post-Inf. Pow. (%\Delta BL)']; [chan_labels{1}, ' - ', chan_labels{2}]}, 'FontSize', 16);
%             
%         end
%         
%         % if g == 1
%         % 
%         %     y = ylabel([band_labels{b}, ' Hz'], 'FontSize', 16, 'Rotation', 0);
%         % 
%         %     set(y, 'Units', 'Normalized', 'Position', [-0.2 0.4 0])
%         % 
%         % end
%         
%         set(gca, 'FontSize', 16)
%         
%         box off
%         
%     end
%     
% end

save_as_eps(gcf, 'STR_M1_kmeans_PLV_by_groups')

save_as_pdf(gcf, 'STR_M1_kmeans_PLV_by_groups')

end
