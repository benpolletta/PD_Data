function PD_GC_pre_post_ttest(filename, folder_indices, band_index)

pd_names = {'Pre-Infusion', 'Post-Infusion'}; no_periods = length(pd_names);

% pd_linestyles = {'--', '-'};

direction_labels = {'Striatum -> Motor Cortex'; 'Motor Cortex -> Striatum'};

no_directions = length(direction_labels);

GC = load(make_sliding_window_analysis_name([filename, '_pre_post_band', num2str(band_index)],...
    'mvgc_analysis', {[150 150], [2 2]}, 2));

GC = abs(GC.GC);

sampling_freq = 500;

freqs = sampling_freq*(1:size(GC, 1))/(2*size(GC,1));

no_freqs = length(freqs);

if isempty(folder_indices)
    
    folder_indices = ones(size(GC, 2), 1);
    
    folder_indices([3 9]) = 0; folder_indices = logical(folder_indices);
    
    folder_indices = {'', folder_indices};
    
end

freq_limit = 100;

p_vals = nan(length(freqs), 2, no_directions);
    
parfor f = 1:no_freqs
    
    temp = nan(2,2);
    
    for direction = 1:no_directions
        
        [~, temp(1, direction)] = ttest(GC(f, folder_indices{2}, 1, direction),...
            GC(f, folder_indices{2}, 2, direction), 'tail', 'left');
        
        [~, temp(2, direction)] = ttest(GC(f, folder_indices{2}, 1, direction),...
            GC(f, folder_indices{2}, 2, direction), 'tail', 'right');
        
    end
    
    p_vals(f, :, :) = reshape(temp, 1, 2, 2); % temp;
    
end

mat_name = make_sliding_window_analysis_name([filename, folder_indices{1}, '_pre_post_band', num2str(band_index)],...
    'mvgc_analysis', {[150 150], [2 2]}, 2);

save([mat_name, '_ttest.mat'], 'p_vals')

