function PD_GC_pre_post_plot(filename, folder_indices, band_index, significance)

figure

pd_names = {'Pre-Infusion', 'Post-Infusion'}; no_periods = length(pd_names);

% pd_linestyles = {'--', '-'};

direction_labels = {'Striatum -> Motor Cortex'; 'Motor Cortex -> Striatum'};

no_directions = length(direction_labels);

analysis_name = make_sliding_window_analysis_name([filename, '_pre_post_band', num2str(band_index)],...
    'mvgc_analysis',{[150 150],[2 2]},2);

load(analysis_name)

sampling_freq = 500;

f = sampling_freq*(1:size(GC, 1))/(2*size(GC,1));

if isempty(folder_indices)
    
    folder_indices = ones(size(GC, 2), 1);
    
    folder_indices([3 9]) = 0; folder_indices = logical(folder_indices);
    
    folder_indices = {'', folder_indices};
    
end

stats_name = [make_sliding_window_analysis_name([filename, folder_indices{1},...
    '_pre_post_band', num2str(band_index)], 'mvgc_analysis', {[150 150],[2 2]}, 2), '_ttest'];

load(stats_name)

test = p_vals < significance;

GC_mean(:, :, :, :) = squeeze(nanmean(GC(:, folder_indices{2}, :, :), 2));
GC_se(:, :, :, :) = squeeze(nanstd(GC(:, folder_indices{2}, :, :), [], 2)/sqrt(sum(folder_indices{2})));

GC_pre_from_post = squeeze(diff(GC(:, folder_indices{2}, :, :), [], 3));

freq_limit = 100;

freq_indices = f <= 100;

for direction = 1:no_directions
    
    h = subplot(2, no_directions, direction);
    
    boundedline(f(freq_indices), GC_mean(freq_indices, :, direction),...
        prep_for_boundedline(GC_se(freq_indices, :, direction))); 
        % squeeze(nanmean(GC_pre_from_post(freq_indices, :, direction), 2)),... 
        % prep_for_boundedline(squeeze(nanstd(GC_pre_from_post(freq_indices, :, direction), [], 2)))); % 
    
    hold on
            
    add_stars(h, f(freq_indices), test(freq_indices, :, direction), [1 0], [1 0 0; 1 .5 0])
    
    set(gca, 'FontSize', 14)
    
    title(direction_labels{direction})
    
    legend(pd_names)
    
    ylabel({'Partial Directed Coherence'; 'Mean \pm S.E.'})
    
    h = subplot(2, no_directions, no_directions + direction);
        
    set(h, 'NextPlot', 'add', 'ColorOrder', distinguishable_colors(size(GC, no_periods)), 'FontSize', 14)
    
    % for pd = 1:2
    
        plot(f(freq_indices), GC_pre_from_post(freq_indices, :, direction)) % , pd_linestyles{pd})
        
    % end
    
    hold on
    
    ylabel({'Partial Directed Coherence'; 'Post-Infusion - Pre-Infusion'})
    
    xlabel('Freq. (Hz)')
    
end

fig_name = make_sliding_window_analysis_name([filename, folder_indices{1},...
    '_pre_post_band', num2str(band_index), '_p', num2str(significance)],...
    'mvgc_analysis',{[150 150],[2 2]},2);

save_as_pdf(gcf, fig_name)

