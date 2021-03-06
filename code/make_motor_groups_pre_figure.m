function make_motor_groups_pre_figure(freq_limit, p_val, test_flag, bands, band_index, period)

pd_labels = {'pre', 'post'}; pd_labels_long = {'Pre-Infusion', 'Post-Infusion'};

pd_index = find(strcmp(pd_labels, period));

load('STR_M1_subjects.mat', 'pd_labels', 'folders')

group_flags = {'M1_not_increased', 'M1_increased'};

group_titles = {'M1+'; 'M1-'};

channel_prefixes = {'STR_w_M1', 'M1'};

chan_labels = {'Striatal', 'M1'};

peak_suffix = '_kmeans';

for b = 1:length(bands);
    
    band_labels{b} = sprintf('%d-%d', bands(b, 1), bands(b,2));
    
end

band_flag = sprintf('%dbands', length(bands));

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

[PLV_mean, PLV_ci] = deal(nan(freq_limit, 2));

for group = 1:2 % Plotting mean and CI.

    load(['STR_M1_1-200Hz_3-21cycles_', band_flag, peak_suffix, '_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_PLV_', group_flags{group} '_Coh_sec_data_for_plot.mat'])
    
    PLV_mean(:, group) = All_mean_mean(1:freq_limit, pd_index);
    
    PLV_ci(:, group) = norminv(1 - max(p_val), 0, 1)*All_mean_se(1:freq_limit, pd_index);
    
end

subplot(3, 2, 5) % 3, 3*(group - 1) + 3); %  2, 4 + group); %

boundedline((1:freq_limit)', PLV_mean, prep_for_boundedline(PLV_ci))

axis tight

set(gca, 'FontSize', 10)

xlabel('Freq. (Hz)', 'FontSize', 10)

ylabel({'PLV'}, 'FontSize', 10) % (% \Delta BL)'; ['Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI']}, 'FontSize', 10); % , 'Rotation', 0);

% set(y, 'Units', 'Normalized', 'Position', [-0.2 0.4 0])

%% Plotting spectra.

for ch = 1:2
    
    % Plotting mean and CI.
    
    [Spec_mean, Spec_ci] = deal(nan(freq_limit, 2));
    
    for group = 1:2
        
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_', band_flag, peak_suffix, '_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
        
        Spec_mean(:, group) = All_mean_mean(1:freq_limit, pd_index);
        
        Spec_ci(:, group) = norminv(1 - max(p_val), 0, 1)*All_mean_se(1:freq_limit, pd_index);
        
    end
    
    subplot(3, 2, (ch - 1)*2 + 1);
    
    boundedline((1:freq_limit)', diag((1:freq_limit).^(3/4))*Spec_mean,...
        prep_for_boundedline(diag((1:freq_limit).^(3/4))*Spec_ci))
    
    axis tight
    
    set(gca, 'FontSize', 10)
    
    ylabel({[chan_labels{ch}, ' Power']}, 'FontSize', 10) % (% \Delta BL)']; ['Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI']}, 'FontSize', 10); % , 'Rotation', 0);
    
    % set(y, 'Units', 'Normalized', 'Position', [-0.2 0.4 0])
    
    if ch == 1
        
        legend(group_titles, 'FontSize', 11, 'Location', 'NE')
        
        title([pd_labels_long{pd_index}, ' Measures'], 'FontSize', 11)
        
    end
    
end

%% Computing stats between pre-infusion spectral power.

for ch = 1:2
    
    load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_spectrum_data_for_plot.mat'])
    
    p_vals = nan(freq_limit, 2);
    
    for pd = pd_index
       
        All_mean(:, :, pd) = All_mean(:, :, pd)*diag(1./sum(All_mean(:, :, pd)));
        
    end
    
    for f = 1:freq_limit % Calculating differences between groups.
        
        if strcmp(test_flag, 'ranksum')
            
            p_vals(f, 1) = ranksum(All_mean(f, logical(folder_index{1}), pd_index)', All_mean(f, logical(folder_index{2}), pd_index)', 'tail', 'left');
            
            p_vals(f, 2) = ranksum(All_mean(f, logical(folder_index{1}), pd_index)', All_mean(f, logical(folder_index{2}), pd_index)', 'tail', 'right');
            
        elseif strcmp(test_flag, 'ttest')
            
            [~, p_vals(f, 1)] = ttest2(All_mean(f, logical(folder_index{1}), pd_index)', All_mean(f, logical(folder_index{2}), pd_index)', 'tail', 'left');
            
            [~, p_vals(f, 2)] = ttest2(All_mean(f, logical(folder_index{1}), pd_index)', All_mean(f, logical(folder_index{2}), pd_index)', 'tail', 'right');
            
        end
        
    end
    
    [test, colors, ~] = test_p_vals(p_vals, p_val, [0 1 1; 1 0 1]); % Getting test matrix, gradient of colors if p_val is a vector.
    
    % test = p_vals < p_val;
    
    mean_mat = nan(freq_limit, 2);
    
    se_mat = nan(freq_limit, 2);
    
    for group = 1:2
        
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
    
        mean_mat(:, group) = All_mean_mean(1:freq_limit, pd_index);
        
        se_mat(:, group) = All_mean_se(1:freq_limit, pd_index);
        
    end
        
    subplot(3, 2, (ch - 1)*2 + 2) % subplot(3, 3, 3*(ch - 1) + 3)
        
    boundedline((1:freq_limit)', -diag((1:freq_limit).^(3/4))*diff(mean_mat, [], 2),...
        prep_for_boundedline(norminv(1 - max(p_val), 0, 1)*diag((1:freq_limit).^(3/4))*sum(sqrt(se_mat).^2, 2)), 'cmap', [0 0 0])
    
    % plot((1:freq_limit)', -diff(mean_mat, [], 2), 'k', 'LineWidth', 1.5)
    
    axis tight % ylim(y_extremes(ch, :))
    
    add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, :, 1)), 0, colors(:, :, 1))
    
    add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, :, 2)), 1, colors(:, :, 2))

    if ~any(test ~= 0), hold on, end
    
    plot((1:freq_limit)', zeros(freq_limit, 1), 'k:')
    
    if ch == 1
        
        title(['(',group_titles{1},') - (', group_titles{2},')'], 'FontSize', 11, 'Rotation', 0)
        
    end
    
    % set(gca, 'FontSize', 16)
    
    box off
    
end

%% Computing stats between pre-infusion PLV.

load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_data_for_plot.mat'])

p_vals = nan(freq_limit, 2);

for f = 1:freq_limit
    
    if strcmp(test_flag, 'ranksum')
        
        p_vals(f, 1) = ranksum(All_mean(f, logical(folder_index{1}), pd_index)', All_mean(f, logical(folder_index{2}), pd_index)', 'tail', 'left');
        
        p_vals(f, 2) = ranksum(All_mean(f, logical(folder_index{1}), pd_index)', All_mean(f, logical(folder_index{2}), pd_index)', 'tail', 'right');
        
    elseif strcmp(test_flag, 'ttest')
        
        [~, p_vals(f, 1)] = ttest2(All_mean(f, logical(folder_index{1}), pd_index)', All_mean(f, logical(folder_index{2}), pd_index)', 'tail', 'left');
        
        [~, p_vals(f, 2)] = ttest2(All_mean(f, logical(folder_index{1}), pd_index)', All_mean(f, logical(folder_index{2}), pd_index)', 'tail', 'right');
        
    end
    
end
        
[test, colors, p_tag] = test_p_vals(p_vals, p_val, [0 1 1; 1 0 1]); % Getting test matrix, gradient of colors if p_val is a vector.

% test = p_vals < p_val;

[mean_mat, se_mat] = deal(nan(freq_limit, 2));

for group = 1:2
    
    load(['STR_M1_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV_', group_flags{group}, '_Coh_sec_data_for_plot.mat'])
    
    mean_mat(:, group) = All_mean_mean(1:freq_limit, pd_index);
        
    se_mat(:, group) = All_mean_se(1:freq_limit, pd_index);
    
end

subplot(3, 2, 6)
        
boundedline((1:freq_limit)', -diff(mean_mat, [], 2),...
    prep_for_boundedline(norminv(1 - max(p_val), 0, 1)*sqrt(sum(se_mat.^2, 2))), 'cmap', [0 0 0])

% plot((1:freq_limit)', -diff(mean_mat(1:freq_limit, :), [], 2), 'k', 'LineWidth', 1.5)

axis tight % ylim(y_extremes(3, :))

add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, :, 1)), 0, colors(:, :, 1))

add_stars_one_line(gca, (1:freq_limit)', logical(test(1:freq_limit, :, 2)), 1, colors(:, :, 2))

if ~any(test ~= 0), hold on, end 

plot((1:freq_limit)', zeros(freq_limit, 1), 'k:')

% set(gca, 'FontSize', 16)

box off

xlabel('Freq. (Hz)', 'FontSize', 11)

%% Saving.

% set(gcf, 'PaperOrientation', 'landscape', 'Units', 'centimeters', 'Position', [0 0 9.1 9.1], 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 9.1 9.1])
% 
% print(gcf, '-painters', '-dpdf', '-r600', [sprintf('STR_M1_kmeans_by_motor_groups_pre_%s_f%d', band_flag, freq_limit), p_tag, '_', test_flag, '.pdf'])
% 
% print(gcf, '-painters', '-depsc', '-r600', [sprintf('STR_M1_kmeans_by_motor_groups_pre_%s_f%d', band_flag, freq_limit), p_tag, '_', test_flag, '.eps'])
% 
% saveas(gcf, [sprintf('STR_M1_kmeans_by_motor_groups_pre_%s_f%d', band_flag, freq_limit), p_tag, '_', test_flag, '.fig'])

save_as_pdf(gcf, [sprintf('STR_M1_kmeans_by_motor_groups_%s_%s_f%d', period, band_flag, freq_limit), p_tag, '_', test_flag])

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
