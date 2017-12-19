function make_supplemental_figure

p_val = .05; freq_limit = 100;

peak_suffix = '_kmeans_win_420_1020';

short_chan_labels = {'M1'};

band_number_flag = {'7', '8', '8', '7', '7', '7'};

band_labels = {'1-4', '5-8', '9-14', '15-30', '40-100', '120-180'}; % , '0-200'};

no_bands = length(band_labels);

load('STR_M1_subjects.mat')

figure

colorspec = [0 0 1; 0 .5 0];

load('missing_2')

%% Plotting period boxplots.

% % This code is in case the times got computed from the
% % STR_M1_subjects.mat file in the same way as all the other
% % subjects.mat files.
% [bp_max_start, bp_max_end] = deal(nan(length(folders), 2, no_bands, 2));
%     
% band_number_labels = {'7', '8'};
% 
% for b = 1:2
%     
%     load(['STR_M1_1-200Hz_3-21cycles_', band_number_labels{b}, 'bands', peak_suffix, '_pct_BP_high_2.5_min_secs_by_STR.mat']) % .mat']) % 
%     
%     basetimes_mat = repmat(basetimes', [1, size(All_bp_max_start, 4), no_bands])/60;
%     
%     for s_id = 1:2
%     
%         subject_index = striatal_id == s_id;
%     
%         bp_max_start(subject_index, :, :, b) = permute(All_bp_max_start(subject_index, s_id, 1:no_bands, :), [1 4 3 2])/(500*60) - basetimes_mat(subject_index, :, :);
%     
%         bp_max_end(subject_index, :, :, b) = permute(All_bp_max_end(subject_index, s_id, 1:no_bands, :), [1 4 3 2])/(500*60) - basetimes_mat(subject_index, :, :);
%     
%     end
%     
% end
    
% This code is in case the times were collected using
% collect_striatal_w_motor_starts_ends.
[bp_max_start, bp_max_end] = deal(nan(length(folders), 2, no_bands, 2));

band_number_labels = {'7', '8'};

for b = 1:2
    
    load(['STR_M1_1-200Hz_3-21cycles_', band_number_labels{b}, 'bands', peak_suffix, '_pct_BP_high_2.5_min_secs_by_STR.mat']) % .mat']) % 
    
    basetimes_mat = repmat(basetimes', [1, size(All_bp_max_start, 4), no_bands])/60;
    
    bp_max_start(:, :, :, b) = permute(All_bp_max_start(:, 1, 1:no_bands, :), [1 4 3 2])/(500*60) - basetimes_mat(:, :, :);
    
    bp_max_end(:, :, :, b) = permute(All_bp_max_end(:, 1, 1:no_bands, :), [1 4 3 2])/(500*60) - basetimes_mat(:, :, :);
    
end

folder_chi = ones(size(folders));

folder_cell = missing_2{2};

for fo = 1:length(folder_cell)
    
    folder_chi(strcmp(folders, folder_cell{fo})) = 0;
    
end

M1_beta = ones(size(folders));

load('M1_groups')

for fo = 1:length(M1_groups{2})
   
    M1_beta(strcmp(folders, M1_groups{2}{fo})) = 0;
    
end

subj_index = find(folder_chi);

mean_bp_max_start = mean(bp_max_start(subj_index, :, :));

std_bp_max_start = std(bp_max_start(subj_index, :, :));

no_measures_plotted = 6;

for b = 1:no_bands

    fprintf('Mean Start of Max. High %s Hz Density = %f\n', band_labels{b}, mean_bp_max_start(:, 2, b))

    fprintf('St. Dev. Start of Max. High %s Hz Density = %f\n', band_labels{b}, std_bp_max_start(:, 2, b))
   
    subplot(no_bands, no_measures_plotted, (b - 1)*no_measures_plotted + 1)
    
    folder_index = 1;
    
    if b == 2 || b == 3, band_number_index = 2; else band_number_index = 1; end
    
    for fo = subj_index % 1:length(folders)
        
        plot([-10 30], folder_index*[1 1], 'k')
        
        hold on
        
        for pd = 1:2
            
            if M1_beta(fo)
                
                plot([bp_max_start(fo, pd, b, band_number_index); bp_max_end(fo, pd, b, band_number_index)], folder_index*[1 1],... % '-d',...
                    'LineWidth', 4.5, 'Color', colorspec(pd, :))
                % 'MarkerFaceColor', colorspec(pd, :), 'MarkerEdgeColor', colorspec(pd, :),...
                
            else
                
                plot([bp_max_start(fo, pd, b, band_number_index); bp_max_end(fo, pd, b, band_number_index)], folder_index*[1 1],... % '-d',...
                    'LineWidth', 2.5, 'Color', colorspec(pd, :)) % , 'LineStyle', '--')
                % 'MarkerFaceColor', colorspec(pd, :), 'MarkerEdgeColor', colorspec(pd, :),...
                
            end
            
        end
        
        folder_index = folder_index + 1;
        
    end
    
    plot([0; 0], [1; length(folders)], 'k')
    
    xlim([-10 30]), ylim([.5 (length(subj_index) +.5)])
    
    box off
    
    if b == 1
        
        title({'Pd. of Highest';'Striatal BPD'}, 'FontSize', 10)
        
    elseif b == no_bands
        
        xlabel({'Time'; '(m, Rel. Infusion)'}, 'FontSize', 10)
        
    end
    
    y = ylabel({band_labels{b}; 'Hz'}, 'FontSize', 10, 'Rotation', 45);
    
    % set(y, 'Units', 'Normalized', 'Position', [-0.35 0.4 0])
    
end

%% Plotting spectra.

chan_prefixes = {'M1'}; no_chans = length(chan_prefixes);

for b = 1:no_bands
    
    for ch = 1:no_chans
        
        load([chan_prefixes{ch}, '_1-200Hz_3-21cycles_', band_number_flag{b}, 'bands', peak_suffix, '_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_missing_2_ch1_data_for_plot.mat'])
        
        All_mean_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
        
        subplot(no_bands, no_measures_plotted, (b - 1)*no_measures_plotted + 1 + ch)
        
        boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))
        
        axis tight
        
        load([chan_prefixes{ch}, '_1-200Hz_3-21cycles_', band_number_flag{b}, 'bands', peak_suffix, '_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat'])
        
        p_vals = nan(freq_limit, 2);
        
        for f = 1:freq_limit
                
                [~, p_vals(f, 1)] = ttest(All_mean(f, logical(folder_chi), 1)', All_mean(f, logical(folder_chi), 2)', 'tail', 'right');
                
                [~, p_vals(f, 2)] = ttest(All_mean(f, logical(folder_chi), 1)', All_mean(f, logical(folder_chi), 2)', 'tail', 'left');
            
        end
        
        test = p_vals < p_val;
        
        % add_stars_one_line(gca, (1:freq_limit)', logical(test(:, 1)), 0, [1 .5 0])
        
        add_stars_one_line(gca, (1:freq_limit)', logical(test(:, 2)), 1, 'c_order', [1 0 0])
        
        if b == 1
            
            title([short_chan_labels{ch}, ' Spectrum'], 'FontSize', 10)
            
        elseif b == no_bands
        
            xlabel('Freq. (Hz)', 'FontSize', 10)
        
        end
        
    end
    
end

%% Plotting PLV.

for b = 1:no_bands
    
    load(['STR_M1_1-200Hz_3-21cycles_', band_number_flag{b}, 'bands', peak_suffix, '_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV_missing_2_Coh_sec_pct_data_for_plot.mat'])
    
    All_mean_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
    
    subplot(no_bands, no_measures_plotted, (b - 1)*no_measures_plotted + 3)
    
    boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))
    
    axis tight
    
    load(['STR_M1_1-200Hz_3-21cycles_', band_number_flag{b}, 'bands', peak_suffix, '_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV_data_for_plot.mat'])
    
    p_vals = nan(freq_limit, 2);
    
    for f = 1:freq_limit
        
        [~, p_vals(f, 1)] = ttest(All_mean(f, logical(folder_chi), 1)', All_mean(f, logical(folder_chi), 2)', 'tail', 'right');
        
        [~, p_vals(f, 2)] = ttest(All_mean(f, logical(folder_chi), 1)', All_mean(f, logical(folder_chi), 2)', 'tail', 'left');
        
    end
    
    test = p_vals < p_val;
    
    % add_stars_one_line(gca, (1:freq_limit)', logical(test(:, 1)), 0, [1 .5 0])
    
    add_stars_one_line(gca, (1:freq_limit)', logical(test(:, 2)), 1, 'c_order', [1 0 0])
    
    if b == 1
        
        title('PLV', 'FontSize', 10)
        
    elseif b == no_bands
        
        xlabel('Freq. (Hz)', 'FontSize', 10)
        
    end
    
end

%% Plotting Cross-Spectrum.
    
pd_names = {'pre', 'post'}; no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods
    
    pd_label = [pd_label, '_', pd_names{period}];
    
end

if nargin < 6, norm_struct = []; end
if isempty(norm_struct), norm_struct = struct('who', 'baseline', 'how', ''); end

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how'), if ~isempty(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how]; end, end

channel_labels = {'Striatum', 'Motor Ctx.'};

for b = 1:no_bands
    
    analysis = 'arxspec'; window_length = 10; band_index = b;
    
    initialize_PD_sliding_window_pre_post
    
    window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);
    
    band_name = [make_sliding_window_analysis_name([filename, pd_label,...
        '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.bands), 7, 'back')], function_name,...
        window_time_cell, 2, varargin{:}), norm_label];
    
    load([band_name, '_recordingspacked.mat'])
    
    mytitle = 'Cross Spectrum';
    
    SW_cross = SW_RecordingsPacked.axissubset('Channel 1', channel_labels{1});
    SW_cross = SW_cross.axissubset('Channel 2', channel_labels{2});
    
    subplot(no_bands, no_measures_plotted, (b - 1)*no_measures_plotted + 4);
    
    xp_comparison_plot_2D(squeeze(SW_cross), @ttest, 2*p_val, [], 0)
    
    legend off
    
    if b == 1
        
        title(mytitle, 'FontSize', 10)
        
    end
    
end

%% Plotting Granger Causality.
    
pd_names = {'pre', 'post'}; no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods
    
    pd_label = [pd_label, '_', pd_names{period}];
    
end

if nargin < 6, norm_struct = []; end
if isempty(norm_struct), norm_struct = struct('who', 'baseline', 'how', ''); end

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how'), if ~isempty(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how]; end, end

channel_labels = {'Striatum', 'Motor Ctx.'};

for b = 1:no_bands
    
    analysis = 'granger'; window_length = 10; band_index = b;
    
    initialize_PD_sliding_window_pre_post
    
    window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

    band_name = [make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.bands), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}), norm_label];
            
    load([band_name, '_recordingspacked.mat'])
    
    channel_loc = [2 1];
    
    for direction = 1:2
        
        channel_loc = fliplr(channel_loc);
        
        direction_title = sprintf('PDC, %s->%s', chan_labels{channel_loc});
        
        SW_direction = SW_RecordingsPacked.axissubset('Channel From', channel_labels{channel_loc(1)});
        SW_direction = SW_direction.axissubset('Channel To', channel_labels{channel_loc(2)});
    
        subplot(no_bands, no_measures_plotted, (b - 1)*no_measures_plotted + 4 + direction);
        
        xp_comparison_plot_2D(squeeze(SW_direction), @ttest, 2*p_val, [], 1)
        
        legend off
        
        if b == 1
        
            title(direction_title, 'FontSize', 10)
            
        end
    
    end
        
end

set(gcf, 'PaperOrientation', 'landscape', 'Units', 'centimeters', 'Position', [0 0 27.3 9.1*1.5], 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 27.3 9.1*1.5])

saveas(gcf, 'supplementary_figure.fig')

print(gcf, '-painters', '-dpdf', '-r600', 'supplementary_figure.pdf')

print(gcf, '-painters', '-depsc', '-r600', 'supplementary_figure.eps')
    
end