function make_motor_groups_figure(freq_limit, p_val, test_flag, bands, band_index, norm_struct)

load('STR_M1_subjects.mat', 'pd_labels', 'folders')

group_flags = {'M1_not_increased', 'M1_increased'};

group_titles = {'M1+'; 'M1-'};

channel_prefixes = {'STR_w_M1', 'M1'};

chan_labels = {'Striatal', 'M1'}; channel_labels = {'Striatum', 'Motor Ctx.'};

no_chans = length(chan_labels);

peak_suffix = '_kmeans';

no_pds_plotted = 2;

for b = 1:length(bands)

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

no_measures_plotted = 6;

ax = nan(no_measures_plotted, 2);

%% Plotting PLV by group.

for group = 1:2 % Plotting mean and CI.

    load(['STR_M1_1-200Hz_3-21cycles_', band_flag, peak_suffix, '_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_PLV_', group_flags{group} '_Coh_sec_pct_data_for_plot.mat'])

    PLV_mean = All_mean_mean; PLV_ci = norminv(1 - max(p_val), 0, 1)*All_mean_se;

    ax(no_measures_plotted, group) = subplot(no_measures_plotted, 2, 4 + group); % no_measures_plotted*(group - 1) + no_measures_plotted) %

    boundedline((1:freq_limit)', PLV_mean(1:freq_limit, 1:no_pds_plotted), prep_for_boundedline(PLV_ci(1:freq_limit, 1:no_pds_plotted)))

    axis tight

    % y_lims(group, :) = ylim;

end

linkaxes(ax(no_measures_plotted, :))

% y_extremes(no_measures_plotted, :) = [min(y_lims(:, 1)) max(y_lims(:, 2))];

for group = 1:2 % Plotting stats.

    % Loading data for all individuals.
    load(['STR_M1_1-200Hz_3-21cycles_', band_flag, '_kmeans_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_PLV_data_for_plot.mat'])

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
    % load(['STR_M1_1-200Hz_3-21cycles_', band_flag, peak_suffix, '_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_PLV_', group_flags{group}, '_Coh_sec_pct_data_for_plot.mat'])
    %
    % PLV_mean = All_mean_mean; PLV_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
    %
    % [sig_lower, sig_higher] = find_sig(PLV_mean(:, 1:2), PLV_ci(:, 1:2));

    subplot(no_measures_plotted, 2, 4 + group) % no_measures_plotted*(group - 1) + no_measures_plotted) %

    % ylim(y_extremes(no_measures_plotted, :))

    add_stars(gca, (1:freq_limit)', logical(test(:, :, 2)), 1, [1 0 0])

    % add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 1)), 0, colors(:, :, 1)) % logical(test(:, 1)), 0, [1 .5 0])

    % add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 2)), 1, colors(:, :, 2)) % logical(test(:, 2)), 1, [1 0 0])

    % add_stars_one_line(gca, (1:freq_limit)', logical(sig_lower(1:freq_limit)), 0, [1 .5 0])
    %
    % add_stars_one_line(gca, (1:freq_limit)', logical(sig_higher(1:freq_limit)), 1, [1 0 0])

    hold on

    % plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')

    % plot((1:freq_limit)', zeros(1, freq_limit), ':k')

    % axis tight

    set(gca, 'FontSize', 10)

    % xlabel('Freq. (Hz)', 'FontSize', 10)

    if group == 1

        y = ylabel({'PLV'}, 'FontSize', 10); % (% \Delta BL)'; ['Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI']}, 'FontSize', 10); % , 'Rotation', 0);

        % set(y, 'Units', 'Normalized', 'Position', [-0.2 0.4 0])

    end

end

%% Plotting spectra.

for ch = 1:2

    % Plotting mean and CI.
    for group = 1:2

        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_', band_flag, peak_suffix, '_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_pct_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])

        All_mean_ci = norminv(1 - max(p_val), 0, 1)*All_mean_se;

        ax(ch, group) = subplot(no_measures_plotted, 2, 2*(ch - 1) + group); % no_measures_plotted*(group - 1) + ch) %

        boundedline((1:freq_limit)', All_mean_mean(1:freq_limit, :), prep_for_boundedline(All_mean_ci(1:freq_limit, :)))

        axis tight

        % y_lims(group, :) = ylim;

    end

    linkaxes(ax(ch, :))

    % y_extremes(ch, :) = [min(y_lims(:, 1)) max(y_lims(:, 2))];

    % Plotting stats.
    for group = 1:2

        % Loading data for all individuals.
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_', band_flag, '_kmeans_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat'])

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

        subplot(no_measures_plotted, 2, 2*(ch - 1) + group) % no_measures_plotted*(group - 1) + ch) %

        % ylim(y_extremes(ch, :))

        % add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 1)), 0, colors(:, :, 1))

        % add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 2)), 1, colors(:, :, 2))

        add_stars(gca, (1:freq_limit)', logical(test(:, :, 2)), 1, [1 0 0])

        % % Calculating non-overlap of CIs.
        % load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_', band_flag, peak_suffix, '_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_pct_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
        %
        % All_mean_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
        %
        % [sig_lower, sig_higher] = find_sig(All_mean_mean, All_mean_ci);
        %
        % add_stars_one_line(gca, (1:freq_limit)', logical(sig_lower(1:freq_limit)), 0, [1 .5 0])
        %
        % add_stars_one_line(gca, (1:freq_limit)', logical(sig_higher(1:freq_limit)), 1, [1 0 0])

        % plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')

        set(gca, 'FontSize', 10)

        % xlabel('Freq. (Hz)', 'FontSize', 10)

        if group == 1

            y = ylabel({[chan_labels{ch}, ' Power']}, 'FontSize', 10); % (% \Delta BL)']; ['Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI']}, 'FontSize', 10); % , 'Rotation', 0);

            % set(y, 'Units', 'Normalized', 'Position', [-0.2 0.4 0])

        end

        if ch == 1

            title(group_titles{group}, 'FontSize', 11)

        end

    end

end

%% Plotting Granger with xPlt.

analysis = 'granger'; window_length = 10;

initialize_PD_sliding_window_pre_post

window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

pd_names = {'pre', 'post'}; no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods

    pd_label = [pd_label, '_', pd_names{period}];

end

if nargin < 6, norm_struct = []; end
if isempty(norm_struct), norm_struct = struct('who', 'baseline', 'how', ''); end

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how'), if ~isempty(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how]; end, end

M1_increased_index{2} = M1_increased_index{2} & All_index{2};

M1_not_increased_index{2} = M1_not_increased_index{2} & All_index{2};

groups_plotted = {M1_increased_index, M1_not_increased_index}; % {All_index};

for group = 1:length(groups_plotted)

    group_name = [make_sliding_window_analysis_name([filename, groups_plotted{group}{1}, pd_label,...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', legnth(data_labels_struct.bands), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}), norm_label];

    load([group_name, '_recordingspacked.mat'])

    frequencies = SW_RecordingsPacked.meta.matrix_dim_1.values;
    freq_index = frequencies <= freq_limit;
    SW_RecordingsPacked.data = cellfun(@(x) x(freq_index, :), SW_RecordingsPacked.data, 'UniformOutput', false);
    SW_RecordingsPacked.meta.matrix_dim_1.values = frequencies(freq_index);

    channel_loc = [2 1];

    for direction = 1:2

        channel_loc = fliplr(channel_loc);

        direction_title = sprintf('%s->%s', chan_labels{channel_loc});

        SW_direction = SW_RecordingsPacked.axissubset('Channel From', channel_labels{channel_loc(1)});
        SW_direction = SW_direction.axissubset('Channel To', channel_labels{channel_loc(2)});
    
        ax(4 + direction, group) = subplot(no_measures_plotted, 2, 2*(4 + direction - 1) + group);

        xp_comparison_plot_2D(squeeze(SW_direction), str2func(test_flag), 2*p_val, [], 1)

        legend off

        if group == 1

            ylabel(direction_title, 'FontSize', 10)

        end

    end

end

sync_axes(ax(5, :)), sync_axes(ax(6, :))

%% Plotting Cross Spectrum with xPlt.

analysis = 'arxspec'; window_length = 10;

initialize_PD_sliding_window_pre_post

window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

pd_names = {'pre', 'post'}; no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods

    pd_label = [pd_label, '_', pd_names{period}];

end

if nargin < 6, norm_struct = []; end
if isempty(norm_struct), norm_struct = struct('who', 'baseline', 'how', ''); end

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how'), if ~isempty(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how]; end, end

M1_increased_index{2} = M1_increased_index{2} & All_index{2};

M1_not_increased_index{2} = M1_not_increased_index{2} & All_index{2};

groups_plotted = {M1_increased_index, M1_not_increased_index}; % {All_index};

for group = 1:length(groups_plotted)

    group_name = [make_sliding_window_analysis_name([filename, groups_plotted{group}{1}, pd_label,...
        '_band', num2str(data_labels_struct.band_index), make_label('bands', legnth(data_labels_struct.bands), 7, 'back')], function_name,...
        window_time_cell, 2, varargin{:}), norm_label];

    load([group_name, '_recordingspacked.mat'])

    frequencies = SW_RecordingsPacked.meta.matrix_dim_1.values;
    freq_index = frequencies <= freq_limit;
    SW_RecordingsPacked.data = cellfun(@(x) x(freq_index, :), SW_RecordingsPacked.data, 'UniformOutput', false);
    SW_RecordingsPacked.meta.matrix_dim_1.values = frequencies(freq_index);

    SW_cross = SW_RecordingsPacked.axissubset('Channel 1', channel_labels{1});
    SW_cross = SW_cross.axissubset('Channel 2', channel_labels{2});

    ax(4, group) = subplot(no_measures_plotted, 2, 2*3 + group);

    xp_comparison_plot_2D(squeeze(SW_cross), str2func(test_flag), 2*p_val, [], 1)

    legend off

    if group == 1

        ylabel('Cross-Spectrum', 'FontSize', 10)

    end

end

sync_axes(ax(4, :))

%% Saving figure.

name = sprintf('STR_M1_kmeans_by_motor_groups_%s_f%d_%sHz_%s_%s', band_flag, freq_limit, band_labels{band_index}, p_tag, test_flag);

set(gcf, 'PaperOrientation', 'landscape', 'Units', 'centimeters', 'Position', [0 0 9.1 18.2], 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 9.1 18.2])

print(gcf, '-painters', '-dpdf', '-r600', [name, '.pdf'])

print(gcf, '-painters', '-depsc', '-r600', [name, '.eps'])

saveas(gcf, [name, '.fig'])

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
