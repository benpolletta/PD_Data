function make_PLV_figure(freq_limit, p_val, test_flag)

load('STR_M1_subjects.mat', 'chan_labels', 'pd_labels', 'folders')

group_flags = {'M1_not_increased', 'M1_increased'};

group_titles = {'M1+'; 'M1-'};

channel_prefixes = {'STR_w_M1', 'M1'};

no_chans = length(chan_labels);

peak_suffix = '_kmeans';

no_pds_plotted = 2;

figure

%% Calculating indicator function of which individuals are included in each group.

load('M1_groups.mat')

for group = 1:2
    
    no_excluded = length(M1_groups{3 - group}); % M1_groups has M1_increased as first entry, M1_not_increased as second entry.
    
    folder_index{group} = ones(1, length(folders));
    
    for e = 1:no_excluded
    
        folder_index{group} = folder_index{group} - strcmp(folders, M1_groups{3 - group}{e});
        
    end
    
end

%% Plotting PLV by group.

for group = 1:2 % Plotting mean and CI.

    load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_', group_flags{group} '_Coh_sec_pct_data_for_plot.mat'])
    
    PLV_mean = All_mean_mean; PLV_ci = norminv(1 - max(p_val), 0, 1)*All_mean_se;
    
    subplot(3, 3, 3*(group - 1) + 3) % 6 + group)
    
    boundedline((1:freq_limit)', PLV_mean(1:freq_limit, 1:no_pds_plotted), prep_for_boundedline(PLV_ci(1:freq_limit, 1:no_pds_plotted)))
    
    axis tight
    
    y_lims(group, :) = ylim;
    
end

y_extremes(3, :) = [min(y_lims(:, 1)) max(y_lims(:, 2))];

for group = 1:2 % Plotting stats.
    
    % Loading data for all individuals.
    load('STR_M1_1-200Hz_3-21cycles_7bands_kmeans_pct_15-30Hz_high_2.5_min_secs_PLV_data_for_plot.mat')
    
    % Calculating difference between pre-infusion and post-infusion.
    p_vals = nan(freq_limit, 2);
    
    for f = 1:freq_limit
        
        if strcmp(test_flag, 'ranksum')
            
            p_vals(f, 1) = ranksum(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'right');
            
            p_vals(f, 2) = ranksum(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'left');
            
        elseif strcmp(test_flag, 'ttest')
            
            [~, p_vals(f, 1)] = ttest(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'right');
            
            [~, p_vals(f, 2)] = ttest(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'left');
            
        end
        
    end
    
    [test, colors, p_tag] = test_p_vals(p_vals, p_val, [1 .5 0; 1 0 0]); % Getting test matrix, gradient of colors if p_val is a vector.
    
    % test = p_vals < p_val;
    
    % % Calculating overlap of CIs.
    % load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_', group_flags{group}, '_Coh_sec_pct_data_for_plot.mat'])
    % 
    % PLV_mean = All_mean_mean; PLV_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
    % 
    % [sig_lower, sig_higher] = find_sig(PLV_mean(:, 1:2), PLV_ci(:, 1:2));
    
    subplot(3, 3, 3*(group - 1) + 3) % 6 + group)
    
    ylim(y_extremes(3, :))
    
    add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 1)), 0, colors(:, :, 1)) % logical(test(:, 1)), 0, [1 .5 0])
    
    add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 2)), 1, colors(:, :, 2)) % logical(test(:, 2)), 1, [1 0 0])
    
    % add_stars_one_line(gca, (1:freq_limit)', logical(sig_lower(1:freq_limit)), 0, [1 .5 0])
    % 
    % add_stars_one_line(gca, (1:freq_limit)', logical(sig_higher(1:freq_limit)), 1, [1 0 0])
    
    hold on
    
    % plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')
    
    % plot((1:freq_limit)', zeros(1, freq_limit), ':k')
    
    % axis tight
    
    set(gca, 'FontSize', 16)
    
    xlabel('Freq. (Hz)', 'FontSize', 16)
    
    if group == 1
    
        title({'Phase-Locking Magnitude'; ['(%\Delta BL, Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI)']}, 'FontSize', 16)
    
    end
    
end

%% Plotting spectra.

for ch = 1:2
    
    % Plotting mean and CI.
    for group = 1:2
        
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
        
        All_mean_ci = norminv(1 - max(p_val), 0, 1)*All_mean_se;
        
        subplot(3, 3, 3*(group - 1) + ch) % 3*(ch - 1) + group)
        
        boundedline((1:freq_limit)', All_mean_mean(1:freq_limit, :), prep_for_boundedline(All_mean_ci(1:freq_limit, :)))
        
        axis tight
        
        y_lims(group, :) = ylim;
        
    end

    y_extremes(ch, :) = [min(y_lims(:, 1)) max(y_lims(:, 2))];
    
    % Plotting stats.
    for group = 1:2
        
        % Loading data for all individuals.
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands_kmeans_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat'])
        
        p_vals = nan(freq_limit, 2);
        
        for f = 1:freq_limit % Calculating differences between pre- and post-infusion.
            
            if strcmp(test_flag, 'ranksum')
                
                p_vals(f, 1) = ranksum(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'right');
                
                p_vals(f, 2) = ranksum(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'left');
                
            elseif strcmp(test_flag, 'ttest')
                
                [~, p_vals(f, 1)] = ttest(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'right');
                
                [~, p_vals(f, 2)] = ttest(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'left');
                
            end
            
        end
        
        [test, colors, p_tag] = test_p_vals(p_vals, p_val, [1 .5 0; 1 0 0]); % Getting test matrix, gradient of colors if p_val is a vector.
        
        % test = p_vals < p_val;
        
        subplot(3, 3, 3*(group - 1) + ch) % 3*(ch - 1) + group)
        
        ylim(y_extremes(ch, :))
        
        add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 1)), 0, colors(:, :, 1))
        
        add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 2)), 1, colors(:, :, 2))
        
        % % Calculating non-overlap of CIs.
        % load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
        % 
        % All_mean_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
        % 
        % [sig_lower, sig_higher] = find_sig(All_mean_mean, All_mean_ci);
        % 
        % add_stars_one_line(gca, (1:freq_limit)', logical(sig_lower(1:freq_limit)), 0, [1 .5 0])
        % 
        % add_stars_one_line(gca, (1:freq_limit)', logical(sig_higher(1:freq_limit)), 1, [1 0 0])
        
        % plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')
        
        set(gca, 'FontSize', 16)
        
        % xlabel('Freq. (Hz)', 'FontSize', 16)
        
        if group == 1
            
            title({[chan_labels{ch}, ' Power']; ['(% \Delta BL, Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI)']}, 'FontSize', 16)
            
        end
        
        if ch == 1
            
            y = ylabel(group_titles{group}, 'FontSize', 20, 'Rotation', 0); % title(chan_labels{1}, 'FontSize', 20)
            
            set(y, 'Units', 'Normalized', 'Position', [-0.35 0.4 0])
            
        end
        
    end
    
end

%% Computing stats between post-infusion increases in spectral power.

for ch = 1:2
    
    load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat'])
    
    p_vals = nan(freq_limit, 2);
    
    for f = 1:freq_limit % Calculating differences between groups.
        
        if strcmp(test_flag, 'ranksum')
            
            p_vals(f, 1) = ranksum(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'left');
            
            p_vals(f, 2) = ranksum(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'right');
            
        elseif strcmp(test_flag, 'ttest')
            
            [~, p_vals(f, 1)] = ttest2(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'left');
            
            [~, p_vals(f, 2)] = ttest2(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'right');
            
        end
        
    end
    
    [test, colors, p_tag] = test_p_vals(p_vals, p_val, [0 1 1; 1 0 1]); % Getting test matrix, gradient of colors if p_val is a vector.
    
    % test = p_vals < p_val;
    
    mean_mat = nan(freq_limit, 2);
    
    se_mat = nan(freq_limit, 2);
    
    for group = 1:2
        
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
    
        mean_mat(:, group) = All_mean_mean(1:freq_limit, 2);
        
        se_mat(:, group) = All_mean_se(1:freq_limit, 2);
        
    end
        
    subplot(3, 3, 3*2 + ch) % subplot(3, 3, 3*(ch - 1) + 3)
        
    boundedline((1:freq_limit)', -diff(mean_mat, [], 2), prep_for_boundedline(norminv(1 - max(p_val), 0, 1)*sqrt(sum(se_mat.^2, 2))), 'cmap', [0 0 0])
    
    % plot((1:freq_limit)', -diff(mean_mat, [], 2), 'k', 'LineWidth', 1.5)
    
    axis tight % ylim(y_extremes(ch, :))
    
    add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, :, 1)), 0, colors(:, :, 1))
    
    add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, :, 2)), 1, colors(:, :, 2))

    if ~any(test ~= 0), hold on, end
    
    plot((1:freq_limit)', zeros(freq_limit, 1), 'k:')
    
    if ch == 1
        
        y = ylabel({group_titles{1}; [' - ', group_titles{2}]}, 'FontSize', 20, 'Rotation', 0);
    
        set(y, 'Units', 'Normalized', 'Position', [-0.4 0.4 0])
        
    end
    
    set(gca, 'FontSize', 16)
    
    box off
    
end

%% Computing stats between post-infusion increases in PLV.

load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_Coh_sec_pct_data_for_plot.mat'])

p_vals = nan(freq_limit, 2);

for f = 1:freq_limit
    
    if strcmp(test_flag, 'ranksum')
        
        p_vals(f, 1) = ranksum(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'left');
        
        p_vals(f, 2) = ranksum(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'right');
        
    elseif strcmp(test_flag, 'ttest')
        
        [~, p_vals(f, 1)] = ttest2(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'left');
        
        [~, p_vals(f, 2)] = ttest2(All_mean(f, logical(folder_index{1}), 2)', All_mean(f, logical(folder_index{2}), 2)', 'tail', 'right');
        
    end
    
end
        
[test, colors, p_tag] = test_p_vals(p_vals, p_val, [0 1 1; 1 0 1]); % Getting test matrix, gradient of colors if p_val is a vector.

% test = p_vals < p_val;

[mean_mat, se_mat] = deal(nan(freq_limit, 2));

for group = 1:2
    
    load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_', group_flags{group}, '_Coh_sec_pct_data_for_plot.mat'])
    
    mean_mat(:, group) = All_mean_mean(1:freq_limit, 2);
        
    se_mat(:, group) = All_mean_se(1:freq_limit, 2);
    
end

subplot(3, 3, 3*2 + 3)
        
boundedline((1:freq_limit)', -diff(mean_mat, [], 2), prep_for_boundedline(norminv(1 - max(p_val), 0, 1)*sqrt(sum(se_mat.^2, 2))), 'cmap', [0 0 0])

% plot((1:freq_limit)', -diff(mean_mat(1:freq_limit, :), [], 2), 'k', 'LineWidth', 1.5)

axis tight % ylim(y_extremes(3, :))

add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, :, 1)), 0, colors(:, :, 1))

add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, :, 2)), 1, colors(:, :, 2))

if ~any(test ~= 0), hold on, end 

plot((1:freq_limit)', zeros(freq_limit, 1), 'k:')

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
%         add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, 1)), 0, [1 1 0])
%         
%         add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, 2)), 1, [.5 0 .5])
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

save_as_eps(gcf, [sprintf('STR_M1_kmeans_PLV_by_groups_f%d', freq_limit), p_tag, '_', test_flag])

save_as_pdf(gcf, [sprintf('STR_M1_kmeans_PLV_by_groups_f%d', freq_limit), p_tag, '_', test_flag])

end

function [test, colors, p_tag] = test_p_vals(p_vals, p_val, colors_in)

p_val = sort(p_val);

no_ps = length(p_val);

no_tests = size(p_vals, 2);

test = nan([size(p_vals, 1), no_ps, no_tests]);

colors = nan(no_ps, 3, no_tests);
    
for t = 1:no_tests
    
    test(:, 1, t) = p_vals(:, t) < p_val(1);
    
    for p = 2:no_ps
    
        test(:, p, t) = p_vals(:, t) >= p_val(p - 1) & p_vals(:, t) < p_val(p);
    
    end

    colors(:, :, t) = flipud(color_gradient(no_ps, .75*colors_in(t, :), colors_in(t, :)));
    
end

if isscalar(p_val)
    
    p_tag = sprintf('_p%g', p_val);
    
else
    
    p_tag = sprintf('_p%gto%g', p_val(1), p_val(end));

end

end
