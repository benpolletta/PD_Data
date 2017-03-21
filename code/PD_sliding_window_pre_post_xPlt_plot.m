function PD_sliding_window_pre_post_xPlt_plot(function_name, sliding_window_cell, data_labels_struct, filename, significance, varargin)
    
% Loads sliding window analysis on carbachol data at times of highest
% striatal beta band density, pre- and post-infusion.
% INPUTS:
%     function_handle (function handle): analysis to be performed on each window.
%     sliding_window_cell (2 x 1 cell of 2 x 1 arrays): cell
%       containing sliding window length and step length (in indices) for the
%       two dimensions of the carbachol data.
%     subjects_mat_cell (n x 1 cell of strings): cell containing names of
%       *subjects.mat files to be (batch) analyzed.
%     data_labels_struct (structure): can be initialized by
%       init_data_labels.m, contains fields: BP_suffix, peak_suffix,
%       data_suffix, epoch_secs, pd_handle, data_field, band_index,
%       sampling_freq.     
%     filename (string): name of collected analysis (e.g., 'STR_w_M1').

%% Loading data & putting into xPlt.
function_name = get_fname(function_name);
    
window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

pd_names = {'pre', 'post'}; no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods
    
    pd_label = [pd_label, '_', pd_names{period}];
    
end

load([make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index)], function_name,...
    window_time_cell, 2, varargin{:}), '_xPlt.mat'])

SW_xPlt = abs(SW_xPlt);

SW_xPlt.getaxisinfo

frequencies = SW_xPlt.meta.matrix_dim_1.values;

if sliding_window_cell{1}(1) > 1000
    
    ds_factor = sliding_window_cell{1}(1)/1000;
    
    ds_factors = factor(ds_factor);
    
    for i = 1:length(ds_factors)

        SW_xPlt.data = cellfun(@(x) decimate(x, ds_factors(i)), SW_xPlt.data, 'UniformOutput', 0);
        
        frequencies = decimate(frequencies, ds_factors(i));
        
    end
    
    % SW_xPlt.data = cellfun(@(x) decimate(x, 10^mod(log(ds_factor)/log(10), 1)), SW_xPlt.data, 'UniformOutput', 0);
    % 
    % frequencies = decimate(frequencies, 10^mod(log(ds_factor)/log(10), 1));
    
end

freq_limit = 50; freq_indicator = frequencies <= freq_limit;

freq_mult = diag(frequencies(freq_indicator).^(2/3));

SW_xPlt.data = cellfun(@(x) freq_mult*x(freq_indicator), SW_xPlt.data, 'UniformOutput', 0);

SW_xPlt.meta.matrix_dim_1.values = frequencies(freq_indicator);

load('M1_groups')

M1_increased_index{2} = M1_increased_index{2} & All_index{2};

M1_not_increased_index{2} = M1_not_increased_index{2} & All_index{2};

groups_plotted = {All_index, M1_increased_index, M1_not_increased_index}; % {All_index}; 

for group = 1:length(groups_plotted)
    
    group_name = make_sliding_window_analysis_name([filename, groups_plotted{group}{1}, pd_label,...
    '_band', num2str(data_labels_struct.band_index)], function_name,...
    window_time_cell, 2, varargin{:});
    
    SW_group = SW_xPlt.packDim('Recording', 2);
    
    folders = SW_group.meta.matrix_dim_2.values;
    
    SW_group.data = cellfun(@(x) x(:, groups_plotted{group}{2}), SW_group.data, 'UniformOutput', 0);
    
    SW_group = SW_group.unpackDim(2);
    
    SW_group.axis(findAxes(SW_group.axis,'Recording')).values = ... % SW_xPlt = SW_xPlt.importAxisValues(SW_xPlt, 'Recording', {folders(groups_plotted{group}{2})});
        folders(groups_plotted{group}{2});

    if ~isempty(SW_group.axis.findAxes('Window_Dim_1'))
        %% Plotting if there are sliding windows.
        
        close('all')
        
        %% % % False color image of all windows by recording. % % %
        
        SW_WindowsPacked = squeeze(SW_group.packDim('Window_Dim_1', 2));
        
        SW_WindowsPacked.getaxisinfo
        
        function_handles = {@xp_subplot_grid_adaptive,@xp_matrix_imagesc};
        function_arguments = {{{'Recording', 'Period', 'Channel'}},{}};
        dimensions = {1:length(size(SW_WindowsPacked)),0};
        recursivePlot(SW_WindowsPacked,function_handles,dimensions,function_arguments);
        
        save_all_figs([group_name, '_windows'])
        
        close('all')
        
        % waitforbuttonpress
        
        %% % % Pre-post comparison over windows, plotted by recording. % % %
        
        function_handles = {@xp_subplot_grid_adaptive,@xp_compare_2D};
        function_arguments = {{},{@ttest, significance}};
        dimensions = {{'Recording', 'Channel'},{'Period'}};
        recursivePlot(SW_WindowsPacked,function_handles,dimensions,function_arguments);
        
        save_all_figs([group_name, '_windows_compared_by_subject'])
        
        close('all')
        
        % waitforbuttonpress
        
        %% % % Pre-post comparison of mean across windows, over recordings. % % %
        
        SW_WindowMean = SW_WindowsPacked;
        
        SW_WindowMean.data = cellfun(@(x) mean(x, 2), SW_WindowMean.data, 'UniformOutput', 0);
        
        SW_WindowMean = squeeze(SW_WindowMean.packDim('Recording', 2));
        
        SW_WindowMean.getaxisinfo
        
        dimensions = {{'Channel'},{'Period'}};
        recursivePlot(SW_WindowMean,function_handles,dimensions,function_arguments);
        
        save_all_figs([group_name, '_windows_compared'])
        
    else
        %% Plotting if there are no sliding windows.
        
        close('all')
        
        %% % % Plots by recording. % % %
        
        SW_RecordingsPacked = squeeze(SW_group.packDim('Recording', 2));
        
        function_handles = {@xp_subplot_grid_adaptive,@xp_matrix_basicplot};
        function_arguments = {{{'Period', 'Channel'}},{}};
        dimensions = {1:length(size(SW_RecordingsPacked)),0};
        recursivePlot(SW_RecordingsPacked,function_handles,dimensions,function_arguments);
        
        save_all_figs([group_name, '_recordings'])
        
        close('all')
        
        % waitforbuttonpress
        
        %% % % Pre-post comparisons over recordings. % % %
        
        function_handles = {@xp_subplot_grid_adaptive,@xp_compare_2D};
        function_arguments = {{},{@ttest, significance}};
        dimensions = {{'Channel'},{'Period'}};
        recursivePlot(SW_RecordingsPacked,function_handles,dimensions,function_arguments);
        
        save_all_figs([group_name, '_recordings_compared'])
        
    end
    
end

end

