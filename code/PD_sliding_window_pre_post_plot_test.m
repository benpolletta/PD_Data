function PD_sliding_window_pre_post_plot_test(function_name, sliding_window_cell, data_labels_struct, filename, significance, norm_struct, varargin)

%% Preliminaries.
function_name = get_fname(function_name);
    
window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

pd_names = {'pre', 'post'}; no_periods = length(pd_names);

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how')
    
    if isstring(norm_struct.how)
        
        if ~isempty(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how]; end
        
    elseif iscellstr(norm_struct.how)
        
        if any(cellfun(@(x) ~isempty(x), norm_struct.how))
            
            for c = 1:length(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how{c}]; end
            
        end
        
    end
    
end

pd_label = '';

for period = 1:no_periods
    
    pd_label = [pd_label, '_', pd_names{period}];
    
end

channel_labels = {'Striatum', 'Motor Ctx.'};
short_channel_labels = {'Str.', 'M1'};

%% Loading & looping over groups.

load('M1_groups')

M1_increased_index{2} = M1_increased_index{2} & All_index{2};

M1_not_increased_index{2} = M1_not_increased_index{2} & All_index{2};

groups_plotted = {All_index, M1_increased_index, M1_not_increased_index}; % {All_index}; 

for group = 1:length(groups_plotted)
    
    group_name = [make_sliding_window_analysis_name([filename, groups_plotted{group}{1}, pd_label,...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', legnth(data_labels_struct.bands), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}), norm_label];

    load([group_name, '_by_subject_windows_compared_p', num2str(significance), '_test.mat'])
    
    frequencies = test.meta.matrix_dim_1.values;
    
    freq_limit = 50; freq_indicator = frequencies <= freq_limit;
    
    frequencies = frequencies(freq_indicator);
    test.meta.matrix_dim_1.values = frequencies;
    
    test.data = cellfun(@(x) double(x(freq_indicator, 1)), test.data, 'UniformOutput', 0);
    test_recmean = mean_over_axis(test, 'Recording');
    
    function_handles = {@xp_tight_subplot_adaptive, @xp_matrix};
    function_arguments = {{},{}};
    dimensions = {{'Channel'},{'Period'}};
    recursivePlot(test_recmean, function_handles, dimensions, function_arguments)
   
    save_as_pdf(gcf, [group_name, '_test'])
    
    channel_loc = [2 1];
    
    test_mat = nan(length(frequencies), 2);
    
    for direction = 1:2
        
        channel_loc = fliplr(channel_loc);
        
        test_direction = test_recmean.axissubset('Channel From', channel_labels{channel_loc(1)});
        test_direction = test_direction.axissubset('Channel To', channel_labels{channel_loc(2)});
        
        test_mat(:, direction) = test_direction.data{1}*length(test.axis('Recording').values);
        
        direction_labels{direction} = sprintf('PDC, %s->%s', short_channel_labels{channel_loc});
        
    end
    
    figure
    
    plot(frequencies, test_mat, 'LineWidth', 2)
    
    box off
    set(gca, 'FontSize', 16)
    xlabel('Freq. (Hz)')
    ylabel({'Animals Showing Increase'})
    legend(direction_labels)
    
    save_as_pdf(gcf, [group_name, '_test_figure'])
    
end