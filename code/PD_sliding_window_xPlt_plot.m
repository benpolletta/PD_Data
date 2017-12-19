function PD_sliding_window_xPlt_plot(function_name, sliding_window_cell, data_labels_struct, filename, pd_names, significance, norm, varargin)
    
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

%% Loading data.
function_name = get_fname(function_name);
    
window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods
    
    pd_label = [pd_label, '_', pd_names{period}];
    
end

if iscellstr(norm)
    
    norm_label = '';
    
    for n = 1:length(norm)
        
        norm_label = [norm_label, '_', norm{n}];
        
    end
    
elseif ischar(norm)
    
    norm_label = norm;
    
end

load([make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.band_index), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}), '_xPlt.mat'])

if ~strcmp(function_name, 'PAC')
    
    SW_xPlt = xp_abs(SW_xPlt);
    
    frequencies = SW_xPlt.meta.matrix_dim_1.values;
    
else
    
    amp_freqs = SW_xPlt.meta.matrix_dim_1.values;
    
    phase_freqs = SW_xPlt.meta.matrix_dim_2.values;
    
end

if ~isempty(regexp(pd_label, 'shuffle', 'once'))
    
    if isempty(SW_xPlt.findaxis('Shuffles'))
        
        SW_xPlt = SW_xPlt.unpackDim(axes_info_struct.odims + 1);
        
    end
    
    shuffle_axis = SW_xPlt.findaxis('Shuffles');
    
    SW_xPlt.axis(shuffle_axis).name = 'Window_Dim_1';
    
end

%% Normalizing.

for n = 1:length(norm)

    SW_xPlt = PD_sliding_window_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, pd_names, norm{n}, varargin{:});
    
end

size_dim2 = cellfun(@(x) size(x, 2), SW_xPlt.data);
    
if ~any(size_dim2(:) > 1)
    
    %% Downsampling and restricting frequency.
    
    if sliding_window_cell{1}(1) > 1000 && isint(sliding_window_cell{1}(1)/1000)
        
        ds_factor = sliding_window_cell{1}(1)/1000;
        
        ds_factors = factor(ds_factor);
        
        for f = 1:length(ds_factors)
            
            SW_xPlt.data = cellfun(@(x) decimate(x, ds_factors(f)), SW_xPlt.data, 'UniformOutput', 0);
            
            frequencies = decimate(frequencies, ds_factors(f));
            
        end
        
    end
    
    freq_limit = 50; freq_indicator = frequencies <= freq_limit;
    
    frequencies = frequencies(freq_indicator);
    
    SW_xPlt.meta.matrix_dim_1.values = frequencies;
    
    SW_xPlt.data = cellfun(@(x) x(freq_indicator), SW_xPlt.data, 'UniformOutput', 0);
    
end

%% Loading & looping over groups.

load('M1_groups')

M1_increased_index{2} = M1_increased_index{2} & All_index{2};

M1_not_increased_index{2} = M1_not_increased_index{2} & All_index{2};

groups_plotted = {All_index, M1_increased_index, M1_not_increased_index}; % {All_index}; 

for group = 1:length(groups_plotted)
    
    group_name = [make_sliding_window_analysis_name([filename, groups_plotted{group}{1}, pd_label,...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.band_index), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}), norm_label];
  
    recording_pack_dim = SW_xPlt.lastNonSingletonDim + 1;

    SW_group = SW_xPlt.packDim('Recording', recording_pack_dim);
    
    folders = SW_group.meta.(['matrix_dim_', num2str(recording_pack_dim)]).values;
    
    indices = cell(1, recording_pack_dim); indices(:) = {':'}; indices(recording_pack_dim) = groups_plotted{group}(2);
    
    SW_group.data = cellfun(@(x) x(indices{:}), SW_group.data, 'UniformOutput', 0);
    
    SW_group = SW_group.unpackDim(recording_pack_dim);
    
    SW_group.axis(findaxis(SW_group,'Recording')).values = ... % SW_xPlt = SW_xPlt.importAxisValues(SW_xPlt, 'Recording', {folders(groups_plotted{group}{2})});
        folders(groups_plotted{group}{2});
    
    if ~any(size_dim2(:) > 1)
    
        if ~isempty(SW_group.findaxis('Window_Dim_1'))
            %% Plotting if there are sliding windows.
            
            close('all')
            
            SW_WindowMean = mean_over_axis(SW_group, 'Window_Dim_1');
            
            %% % % False color image of all windows by recording. % % %
            
            SW_WindowsPacked = squeeze(SW_group.packDim('Window_Dim_1', 2));
            
            SW_WindowsPacked.getaxisinfo
            
            function_handles = {@xp_tight_subplot_adaptive,@xp_matrix_imagesc};
            function_arguments = {{{'Recording', 'Period', 'Channel'}},{}};
            dimensions = {1:length(size(SW_WindowsPacked)),0};
            recursivePlot(SW_WindowsPacked,function_handles,dimensions,function_arguments);
            
            save_all_figs([group_name, '_windows'])
            
            close('all')
            
            % waitforbuttonpress
            
            %% % % Pre-post comparison over windows, plotted by recording. % % %
            
            function_handles = {@xp_tight_subplot_adaptive,@xp_comparison_plot_2D};
            function_arguments = {{},{@ttest, significance, [], 1}};
            dimensions = {{'Recording', 'Channel'},{'Period'}};
            recursivePlot(SW_WindowsPacked,function_handles,dimensions,function_arguments);
            
            save_all_figs([group_name, '_windows_by_subject_compared_p', num2str(significance)])
            
            close('all')
            
            % waitforbuttonpress
            
            %% % % Pre-post comparison of mean across windows, over recordings. % % %
            
            SW_WindowMean = squeeze(SW_WindowMean.packDim('Recording', 2));
            
            SW_WindowMean = norm_axis_by_value(SW_WindowMean, 'Period', 'pre');
            
            SW_WindowMean.getaxisinfo
            
            dimensions = {{'Channel'},{'Period'}};
            recursivePlot(SW_WindowMean,function_handles,dimensions,function_arguments);
            
            save_all_figs([group_name, '_windows_compared_p', num2str(significance)])
            
        else
            %% Plotting if there are no sliding windows.
            
            close('all')
            
            %% % % Plots by recording. % % %
            
            SW_RecordingsPacked = squeeze(SW_group.packDim('Recording', 2));
            
            function_handles = {@xp_tight_subplot_adaptive,@xp_matrix_basicplot};
            function_arguments = {{{'Period', 'Channel'}},{}};
            dimensions = {1:length(size(SW_RecordingsPacked)),0};
            recursivePlot(SW_RecordingsPacked,function_handles,dimensions,function_arguments);
            
            save_all_figs([group_name, '_recordings'])
            
            close('all')
            
            % waitforbuttonpress
            
            %% % % Pre-post comparisons over recordings. % % %
            
            function_handles = {@xp_tight_subplot_adaptive,@xp_comparison_plot_2D};
            function_arguments = {{},{@ttest, significance, [], 1}};
            dimensions = {{'Channel'},{'Period'}};
            recursivePlot(SW_RecordingsPacked,function_handles,dimensions,function_arguments);
            
            save_all_figs([group_name, '_recordings_compared_p', num2str(significance)])
            
        end
        
    else
        
        if ~isempty(SW_group.findaxis('Window_Dim_1'))
            %% Plotting if there are sliding windows.
            
            close('all')
            
            %% % % Mean over windows, plotted by recording. % % %
            
            SW_WindowMean = mean_over_axis(SW_group, 'Window_Dim_1');
            
            SW_WindowMean.getaxisinfo
            
            function_handles = {@xp_tight_subplot_adaptive,@xp_matrix_imagesc};
            function_arguments = {{[], [], [], 'row'},{1,1}};
            dimensions = {{'Recording', 'Period', 'Channel'},{}};
            recursivePlot(SW_WindowMean,function_handles,dimensions,function_arguments);
            
            save_all_figs([group_name, '_window_mean_by_subject'])
            
            close('all')
            
            %% % % Pre-post comparison over windows, plotted by recording. % % %
            
            SW_WindowsPacked = packDim(SW_group, 'Window_Dim_1');
            
            SW_WindowsPacked.getaxisinfo
            
            function_handles = {@xp_tight_subplot_adaptive,@xp_compare_3D};
            function_arguments = {{},{@ttest,significance,1,'pcolor'}};
            dimensions = {{'Recording', 'Channel'},{'Period'}};
            recursivePlot(SW_WindowsPacked,function_handles,dimensions,function_arguments)
            
            save_all_figs([group_name, '_by_subject_compared_p', num2str(significance)])
            
            close('all')
            
            %% % % Mean across recordings. % % %
            
            SW_RecordingMean = mean_over_axis(SW_WindowMean, 'Recording');
            
            SW_RecordingMean.getaxisinfo
            
            function_handles = {@xp_tight_subplot_adaptive,@xp_matrix_imagesc};
            function_arguments = {{[], [], [], 'row'},{1,1}};
            dimensions = {{'Channel','Period'},{}};
            recursivePlot(SW_RecordingMean,function_handles,dimensions,function_arguments);
            
            save_all_figs([group_name, '_window_mean'])
            
            close('all')
            
            %% % % Pre-post comparison over windows, plotted by recording. % % %
            
            SW_RecordingsPacked = packDim(SW_WindowMean, 'Recording');
            
            SW_RecordingsPacked.getaxisinfo
            
            function_handles = {@xp_tight_subplot_adaptive,@xp_compare_3D};
            function_arguments = {{},{@ttest,significance,1,'pcolor'}};
            dimensions = {{'Channel'},{'Period'}};
            recursivePlot(SW_RecordingsPacked,function_handles,dimensions,function_arguments)
            
            save_all_figs([group_name, '_compared_p', num2str(significance)])
            
            close('all')
            
        else
            %% Plotting if there are no sliding windows.
            
            close('all')
            
            %% % % Plots by recording. % % %
            
            function_handles = {@xp_tight_subplot_adaptive,@xp_matrix_imagesc};
            function_arguments = {{[], [], [], 'row'},{1,1}};
            dimensions = {{'Recording', 'Period', 'Channel'},{}};
            recursivePlot(SW_group,function_handles,dimensions,function_arguments);
            
            save_all_figs([group_name, '_by_subject'])
            
            close('all')
            
            %% % % Pre-post comparisons over recordings. % % %
            
            SW_RecordingMean = mean_over_axis(SW_group, 'Recording');
            
            SW_RecordingMean.getaxisinfo
            
            dimensions = {{'Channel','Period'},{}};
            recursivePlot(SW_RecordingsPacked,function_handles,dimensions,function_arguments);
            
            save_all_figs([group_name, '_subject_mean'])
            
        end
        
    end
    
end

end