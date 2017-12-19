function make_all_band_figures_v2(missing)

group_prefix = 'STR_M1';

short_chan_labels = {'Str.', 'M1'};

band_labels = {'1-4', '4-8', '8-12', '15-30', '40-100', '120-180'}; % , '0-200'};

no_bands = length(band_labels);

load([group_prefix, '_subjects.mat'])

figure

colorspec = [0 0 1; 0 .5 0];

%% Plotting period boxplots.
    
load([group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_pct_BP_high_2.5_min_secs.mat']) % _by_STR.mat'])

basetimes_mat = repmat(basetimes', [1, size(All_bp_max_start, 4), size(All_bp_max_start, 3)])/60;

[bp_max_start, bp_max_end] = deal(nan(length(folders), 2, no_bands + 1));

for s_id = 1:2
    
    s_ind = striatal_id == s_id;
    
    bp_max_start(s_ind, :, :) = permute(All_bp_max_start(s_ind, s_id, :, :), [1 4 3 2])/(500*60) - basetimes_mat(s_ind, :, :);

    bp_max_end(s_ind, :, :) = permute(All_bp_max_end(s_ind, s_id, :, :), [1 4 3 2])/(500*60) - basetimes_mat(s_ind, :, :);
    
end

folder_chi = ones(size(folders));

folder_cell = missing{2};

for fo = 1:length(folder_cell)
    
    folder_chi(strcmp(folders, folder_cell{fo})) = 0;
    
end

subj_index = find(folder_chi);

mean_bp_max_start = mean(bp_max_start(subj_index, :, :));

std_bp_max_start = std(bp_max_start(subj_index, :, :)); 

for b = 1:no_bands

    fprintf('Mean Start of Max. High %s Hz Density = %f\n', band_labels{b}, mean_bp_max_start(:, 2, b))

    fprintf('St. Dev. Start of Max. High %s Hz Density = %f\n', band_labels{b}, std_bp_max_start(:, 2, b))
   
    subplot(no_bands, 4, (b - 1)*4 + 1)
    
    folder_index = 1;
    
    for fo = subj_index % 1:length(folders)
        
        plot([-10 30], folder_index*[1 1], 'k')
        
        hold on
        
        for pd = 1:2
            
            plot([bp_max_start(fo, pd, b); bp_max_end(fo, pd, b)], folder_index*[1 1],... % '-d',...
                'LineWidth', 4, 'Color', colorspec(pd, :))
                % 'MarkerFaceColor', colorspec(pd, :), 'MarkerEdgeColor', colorspec(pd, :),...
                
            
        end
        
        folder_index = folder_index + 1;
        
    end
    
    plot([0; 0], [1; length(folders)], ':k')
    
    xlim([-10 30]), ylim([.5 (length(subj_index) +.5)])
    
    box off
    
    if b == 1
        
        title({'Pd. Densest';'Str. BP'}, 'FontSize', 16)
        
    elseif b == no_bands
        
        xlabel({'Time'; '(m, Rel. Infusion)'}, 'FontSize', 16)
        
    end
    
    y = ylabel([band_labels{b}, ' Hz'], 'FontSize', 16, 'Rotation', 0);
    
    set(y, 'Units', 'Normalized', 'Position', [-0.35 0.4 0])
    
end

no_chans = length(chan_labels);

%% Plotting spectra.

group_prefixes = {'STR_w_M1', 'M1'};

for b = 1:no_bands
    
    for ch = 1:no_chans
        
        load([group_prefixes{ch}, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum', missing{1}, '_ch1_data_for_plot.mat'])
        
        All_mean_ci = norminv(1 - .05, 0, 1)*All_mean_se;
        
        [sig_lower, sig_higher] = find_sig(All_mean_mean, All_mean_ci);
        
        subplot(no_bands, 4, (b - 1)*4 + 1 + ch)
        
        boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))
        
        axis tight
        
        add_stars(gca, (1:200)', logical(sig_lower), 0, [1 .5 0])
        
        add_stars(gca, (1:200)', logical(sig_higher), 1, [1 0 0])
        
        if b == 1
            
            title({[short_chan_labels{ch}, ' Pow. (% \Delta BL)']; 'Mean \pm 95% CI'}, 'FontSize', 16)
            
        elseif b == no_bands
        
            xlabel('Freq. (Hz)', 'FontSize', 16)
        
        end
        
    end
    
end

%% Plotting PLV.

for b = 1:no_bands
    
    load([group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV', missing{1}, '_Coh_sec_pct_data_for_plot.mat'])
    
    % All_mean_ci = norminv(1 - .05, 0, 1)*All_mean_se;
    
    [sig_lower, sig_higher] = find_sig(All_mean_mean, All_mean_ci);
    
    subplot(no_bands, 4, (b - 1)*4 + 4)
    
    boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))
    
    axis tight
    
    add_stars(gca, (1:200)', logical(sig_lower), 0, [1 .5 0])
    
    add_stars(gca, (1:200)', logical(sig_higher), 1, [1 0 0])
    
    if b == 1
        
        title({'Phase-Locking Mag. (% \Delta BL)'; 'Mean \pm 95% CI'}, 'FontSize', 16)
        
    elseif b == no_bands
        
        xlabel('Freq. (Hz)', 'FontSize', 16)
        
    end
    
end

%% Plotting GC.

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

for b = 1:no_bands
        
    band_name = [make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.band_index), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}), norm_label];
            
    load([group_name, '_recordingspacked.mat'])
    
    channel_loc = [2 1];
    
    for direction = 1:2
        
        channel_loc = fliplr(channel_loc);
        
        direction_title = sprintf('%s->%s', chan_labels{channel_loc});
        
        SW_direction = SW_RecordingsPacked.axissubset('Channel From', channel_labels{channel_loc(1)});
        SW_direction = SW_direction.axissubset('Channel To', channel_labels{channel_loc(2)});
    
        ax(4 + direction, group) = subplot(no_measures_plotted, 2, 2*(4 + direction - 1) + group);
        
        xp_compare_2D(squeeze(SW_direction), str2func(test_flag), 2*p_val, [], 1)
        
        legend off
        
        if group == 1
        
            ylabel(direction_title, 'FontSize', 10)
            
        end
    
    end
        
end

save_as_eps(gcf, [group_prefix, '_', missing{1}, '_all_bands_v2'])

save_as_pdf(gcf, [group_prefix, '_', missing{1}, '_all_bands_v2'])

end