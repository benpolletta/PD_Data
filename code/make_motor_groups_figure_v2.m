function make_motor_groups_figure_v2(freq_limits, p_val, test_flag, bands, band_index, norm_struct)

% freq_limits: a matrix, each row containing a 2 vector of limits on frequencies.

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

no_freq_limits = size(freq_limits, 1);

rows = no_measures_plotted;
columns = 2*no_freq_limits;

tight_subplot_handles = tight_subplot(rows, columns);

ax = nan(rows, columns);

for fl = 1:no_freq_limits
    
    frequencies = 1:200;
    
    freq_index = frequencies >= freq_limits(fl, 1) & frequencies <= freq_limits(fl, 2);
    
    %% Plotting PLV by group.
    
    for group = 1:2 % Plotting mean and CI.
        
        load(['STR_M1_1-200Hz_3-21cycles_', band_flag, peak_suffix, '_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_PLV_', group_flags{group} '_Coh_sec_pct_data_for_plot.mat'])
        
        PLV_mean = All_mean_mean; PLV_ci = norminv(1 - max(p_val), 0, 1)*All_mean_se;
        
        row = (group - 1)*no_freq_limits + fl;
        
        ax(3, row) = tight_subplot_handles(2*columns + row);
        axes(ax(3, row))
        
        boundedline(frequencies(freq_index)', PLV_mean(freq_index, 1:no_pds_plotted), prep_for_boundedline(PLV_ci(freq_index, 1:no_pds_plotted)))
        
        axis tight
        
    end
    
end

linkaxes(ax(3, :), 'y')

for fl = 1:no_freq_limits
    
    frequencies = 1:200;
    
    freq_index = frequencies >= freq_limits(fl, 1) & frequencies <= freq_limits(fl, 2);
    
    for group = 1:2 % Plotting stats.
        
        % Loading data for all individuals.
        load(['STR_M1_1-200Hz_3-21cycles_', band_flag, '_kmeans_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_PLV_data_for_plot.mat'])
        
        % Calculating difference between pre-infusion and post-infusion.
        p_vals = nan(sum(freq_index), 2);
        
        for f = frequencies(freq_index)
            
            if strcmp(test_flag, 'ranksum')
                
                p_vals(f, 1) = ranksum(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'right');
                
                p_vals(f, 2) = ranksum(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'left');
                
            elseif strcmp(test_flag, 'ttest')
                
                [~, p_vals(f, 1)] = ttest(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'right');
                
                [~, p_vals(f, 2)] = ttest(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'left');
                
            end
            
        end
        
        [test, colors, p_tag] = test_p_vals(p_vals, p_val, [1 .5 0; 1 0 0]); % Getting test matrix, gradient of colors if p_val is a vector.
        
        row = (group - 1)*no_freq_limits + fl;
        
        axes(tight_subplot_handles(2*columns + row))
        
        add_stars(gca, frequencies(freq_index)', logical(test(:, :, 2)), 1, [1 0 0])
        
        hold on
        
        if ~exist('y_lims', 'var'), y_lims = ylim;
            
        else curr_y_lims = ylim; y_lims = [min(curr_y_lims(1), y_lims(1)) max(curr_y_lims(2), y_lims(2))];
            
        end
        
        plot([8 15 30; 8 15 30], repmat(y_lims', 1, 3), 'k:')
        
        set(gca, 'FontSize', 10, 'XTickLabelMode', 'auto', 'YTickLabelMode', 'auto')
        
        if group == 1 && fl == 1
            
            ylabel({'PLV'}, 'FontSize', 10);
            
        end
        
        if fl == 2
            
            set(gca, 'YTick', [], 'YColor', 'w')
            
        end
        
    end
    
end

clear y_lims
    
%% Plotting spectra.
        
for ch = 1:2
    
    % Plotting mean and CI.
    for fl = 1:no_freq_limits
        
        frequencies = 1:200;
        
        freq_index = frequencies >= freq_limits(fl, 1) & frequencies <= freq_limits(fl, 2);
        
        for group = 1:2
            
            load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_', band_flag, peak_suffix, '_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_pct_spectrum_', group_flags{group}, '_ch1_data_for_plot.mat'])
            
            All_mean_ci = norminv(1 - max(p_val), 0, 1)*All_mean_se;
            
            row = (group - 1)*no_freq_limits + fl;
            
            ax(ch, row) = tight_subplot_handles((ch - 1)*columns + row);
            axes(ax(ch, row))
            
            boundedline(frequencies(freq_index)', All_mean_mean(freq_index, :), prep_for_boundedline(All_mean_ci(freq_index, :)))
            
            axis tight
            
        end
        
    end
    
    linkaxes(ax(ch, :), 'y')
    
    % Plotting stats.
    for fl = 1:no_freq_limits
        
        frequencies = 1:200;
        
        freq_index = frequencies >= freq_limits(fl, 1) & frequencies <= freq_limits(fl, 2);
        
        for group = 1:2
            
            % Loading data for all individuals.
            load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_', band_flag, '_kmeans_pct_', band_labels{band_index}, 'Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat'])
            
            p_vals = nan(sum(freq_index), 2);
            
            for f = frequencies(freq_index) % Calculating differences between pre- and post-infusion.
                
                if strcmp(test_flag, 'ranksum')
                    
                    p_vals(f, 1) = ranksum(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'right');
                    
                    p_vals(f, 2) = ranksum(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'left');
                    
                elseif strcmp(test_flag, 'ttest')
                    
                    [~, p_vals(f, 1)] = ttest(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'right');
                    
                    [~, p_vals(f, 2)] = ttest(All_mean(f, logical(folder_index{group}), 1)', All_mean(f, logical(folder_index{group}), 2)', 'tail', 'left');
                    
                end
                
            end
            
            [test, colors, p_tag] = test_p_vals(p_vals, p_val, [1 .5 0; 1 0 0]);
            
            row = (group - 1)*no_freq_limits + fl;
            
            axes(tight_subplot_handles((ch - 1)*columns + row))
            
            add_stars(gca, frequencies(freq_index)', logical(test(:, :, 2)), 1, [1 0 0])
            
            set(gca, 'FontSize', 10, 'XTickLabelMode', 'auto', 'YTickLabelMode', 'auto')
            
            hold on
            
            if ~exist('y_lims', 'var'), y_lims = ylim;
                
            else curr_y_lims = ylim; y_lims = [min(curr_y_lims(1), y_lims(1)) max(curr_y_lims(2), y_lims(2))];
                
            end
            
            plot([8 15 30; 8 15 30], repmat(y_lims', 1, 3), 'k:')
            
            if group == 1 && fl == 1
                
                y = ylabel({[chan_labels{ch}, ' Power']}, 'FontSize', 10);
                
            end
            
            if fl == 2
                
                set(gca, 'YTick', [], 'YColor', 'w')
                
            end
            
            if ch == 1 && fl == 1
                
                title(group_titles{group}, 'FontSize', 11)
                
            end
            
        end
        
    end
    
    clear y_lims
    
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

for fl = 1:no_freq_limits
    
    for group = 1:length(groups_plotted)
        
        row = (group - 1)*no_freq_limits + fl;
        
        group_name = [make_sliding_window_analysis_name([filename, groups_plotted{group}{1}, pd_label,...
            '_band', num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), norm_label];
        
        load([group_name, '_recordingspacked.mat'])
        
        frequencies = SW_RecordingsPacked.meta.matrix_dim_1.values;
        freq_index = frequencies >= freq_limits(fl, 1) & frequencies <= freq_limits(fl, 2);
        SW_RecordingsPacked.data = cellfun(@(x) x(freq_index, :), SW_RecordingsPacked.data, 'UniformOutput', false);
        SW_RecordingsPacked.meta.matrix_dim_1.values = frequencies(freq_index);
        
        channel_loc = [2 1];
        
        for direction = 1:2
            
            channel_loc = fliplr(channel_loc);
            
            direction_title = sprintf('%s->%s', chan_labels{channel_loc});
            
            SW_direction = SW_RecordingsPacked.axissubset('Channel From', channel_labels{channel_loc(1)});
            SW_direction = SW_direction.axissubset('Channel To', channel_labels{channel_loc(2)});
            
            ax(4 + direction, row) = tight_subplot_handles((4 + direction - 1)*columns + row);
            axes(ax(4 + direction, row))
            
            xp_comparison_plot_2D(squeeze(SW_direction), str2func(test_flag), 2*p_val, [], 1)
            
            legend off
            
            hold on
            
            if ~exist('y_lims', 'var'), y_lims = ylim;
                
            else curr_y_lims = ylim; y_lims = [min(curr_y_lims(1), y_lims(1)) max(curr_y_lims(2), y_lims(2))];
                
            end
            
            plot([8 15 30; 8 15 30], repmat(y_lims', 1, 3), 'k:')
            
            if group == 1 && fl == 1
                
                ylabel(direction_title, 'FontSize', 10)
                
            end
            
            if fl == 2
                
                set(gca, 'YTick', [], 'YColor', 'w')
                
            end
            
        end
        
    end
    
end

clear y_lims

sync_axes(ax(5, :), 'y'), sync_axes(ax(6, :), 'y')
    
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

for fl = 1:no_freq_limits
    
    for group = 1:length(groups_plotted)
        
        row = (group - 1)*no_freq_limits + fl;
        
        group_name = [make_sliding_window_analysis_name([filename, groups_plotted{group}{1}, pd_label,...
            '_band', num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), norm_label];
        
        load([group_name, '_recordingspacked.mat'])
        
        frequencies = SW_RecordingsPacked.meta.matrix_dim_1.values;
        freq_index = frequencies >= freq_limits(fl, 1) & frequencies <= freq_limits(fl, 2);
        SW_RecordingsPacked.data = cellfun(@(x) x(freq_index, :), SW_RecordingsPacked.data, 'UniformOutput', false);
        SW_RecordingsPacked.meta.matrix_dim_1.values = frequencies(freq_index);
        
        SW_cross = SW_RecordingsPacked.axissubset('Channel 1', channel_labels{1});
        SW_cross = SW_cross.axissubset('Channel 2', channel_labels{2});
        
        ax(4, row) = tight_subplot_handles(3*columns + row);
        axes(ax(4, row))
        
        xp_comparison_plot_2D(squeeze(SW_cross), str2func(test_flag), 2*p_val, [], 1)
        
        legend off
        
        hold on
        
        if ~exist('y_lims', 'var'), y_lims = ylim;
            
        else curr_y_lims = ylim; y_lims = [min(curr_y_lims(1), y_lims(1)) max(curr_y_lims(2), y_lims(2))];
            
        end
        
        plot([8 15 30; 8 15 30], repmat(y_lims', 1, 3), 'k:')
        
        if group == 1 && fl == 1
            
            ylabel('Cross-Spectrum', 'FontSize', 10)
            
        end
        
        if fl == 2
            
            set(gca, 'YTick', [], 'YColor', 'w')
            
        end
        
    end
    
end

clear y_lims

sync_axes(ax(4, :), 'y')

%% Saving figure.

name = [sprintf('STR_M1_kmeans_by_motor_groups_%s_allfreqs_%sHz', band_flag, band_labels{band_index}), p_tag, '_', test_flag];

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
