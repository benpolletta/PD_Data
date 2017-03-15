function PD_GC_pre_post_figure(significance, run_stats_flag)

filename = 'STR_w_M1';

band_index = 4;

if isempty(significance), significance = .025; end

figure

pd_names = {'Pre-Infusion', 'Post-Infusion'}; no_periods = length(pd_names);

% pd_linestyles = {'--', '-'};

direction_labels = {'\bf{Striatum -> Motor Cortex}'; '\bf{Motor Cortex -> Striatum}'};

no_directions = length(direction_labels);

analysis_name = make_sliding_window_analysis_name([filename, '_pre_post_band', num2str(band_index)],...
    'mvgc_analysis',{[150 150],[2 2]},2);

load(analysis_name)

sampling_freq = 500;

f = sampling_freq*(1:size(GC, 1))/(2*size(GC,1));

freq_limit = 50;

freq_indices = f <= 50;

freq_mult = diag(f(freq_indices).^(2/3));

folder_index_cell = make_folder_index_cell;

no_folder_indices = length(folder_index_cell);

group_titles = {{'{\bf All Data}';'{\fontsize{10}Partial Directed Coherence}';'{\fontsize{10}Mean \pm SE}'},...
    {'{\bf M1+}';'{\fontsize{10}PDC (Mean \pm SE)}'}, {'{\bf M1-}';'{\fontsize{10}PDC (Mean \pm SE)}'}};

for fi = 1:no_folder_indices
    
    folder_indices = folder_index_cell{fi};
    
    % GC_pre_from_post = squeeze(diff(GC(:, folder_indices{2}, :, :), [], 3));
    
    GC_plot = abs(GC(:, folder_indices{2}, :, :)); % GC_plot = GC_pre_from_post; % 
    
    GC_mean = squeeze(nanmean(GC_plot, 2));
    GC_se = squeeze(nanstd(GC_plot, [], 2)/sqrt(sum(folder_indices{2})));
    GC_ci = norminv(1 - significance)*GC_se;
    
    for direction = 1:no_directions
        
        ax(fi, direction) = subplot(no_folder_indices, no_directions, (fi - 1)*no_directions + direction);
        
        boundedline(f(freq_indices), freq_mult*GC_mean(freq_indices, :, direction),...
            prep_for_boundedline(freq_mult*GC_ci(freq_indices, :, direction)));
        % squeeze(nanmean(GC_pre_from_post(freq_indices, :, direction), 2)),...
        % prep_for_boundedline(squeeze(nanstd(GC_pre_from_post(freq_indices, :, direction), [], 2)))); %
        
        axis tight
        
        hold on
        
    end
    
    sync_axes(ax(fi, :))
    
    if run_stats_flag
    
        PD_GC_pre_post_ttest(filename, folder_indices, band_index)
    
    end
        
    stats_name = [make_sliding_window_analysis_name([filename, folder_indices{1},...
        '_pre_post_band', num2str(band_index)], 'mvgc_analysis', {[150 150],[2 2]}, 2), '_ttest'];
    
    load(stats_name)
    
    test = p_vals < significance;
    
    for direction = 1:no_directions
        
        subplot(no_folder_indices, no_directions, (fi - 1)*no_directions + direction)
        
        add_stars(ax(fi, direction), f(freq_indices), test(freq_indices, :, direction), [1 0], [1 0 0; 1 .5 0])
        
        set(gca, 'FontSize', 14)
        
        if fi == 1
            
            title(direction_labels{direction})
            
        end
        
        if fi == 1 && direction == 1
        
            legend(pd_names)
            
        end
        
        if direction == 1
        
            ylabel(group_titles{fi})
            
        end
        
        if fi == no_folder_indices
        
            xlabel('Freq. (Hz)')
            
        end
    
    end

end

fig_name = make_sliding_window_analysis_name([filename,...
    '_pre_post_band', num2str(band_index), '_p', num2str(significance)],...
    'mvgc_analysis', {[150 150], [2 2]}, 2);

save_as_pdf(gcf, [fig_name, '_altogether'])

end

function folder_index_cell = make_folder_index_cell

folder_indices = ones(14, 1);

folder_indices([3 9]) = 0; folder_indices = logical(folder_indices);

folder_indices = {'', folder_indices};

folder_index_cell{1} = folder_indices;

load('M1_groups')

folder_index_cell{2} = {M1_increased_index{1}, ~M1_not_increased_index{2}};

folder_index_cell{3} = {M1_not_increased_index{1}, ~M1_increased_index{2}};

end

