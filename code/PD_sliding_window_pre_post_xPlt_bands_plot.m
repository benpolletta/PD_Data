function PD_sliding_window_pre_post_xPlt_bands_plot(function_name, sliding_window_cell, data_labels_struct, filename, significance, norm_struct, varargin)

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

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how')
        
    if ~isempty(norm_struct.how)
        
        if ischar(norm_struct.how)
            
            norm_label = [norm_label, '_', norm_struct.how];
            
        elseif iscellstr(norm_struct.how)
            
            if any(cellfun(@(x) ~isempty(x), norm_struct.how))
                
                for c = 1:length(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how{c}]; end
                
            end
            
        end
        
    end
    
end

load([make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index)], function_name,...
    window_time_cell, 2, varargin{:}), '_xPlt.mat'])

SW_xPlt.getaxisinfo

SW_Bands = SW_xPlt;

if strcmp(function_name, 'PAC')

    for i = 1:2,

        frequencies{i} = SW_xPlt.meta.(['matrix_dim_', num2str(i)]).values;

    end

else

    SW_xPlt = xp_abs(SW_xPlt);

    frequencies = SW_xPlt.meta.matrix_dim_1.values;

end

%% Normalizing.

SW_xPlt = PD_sliding_window_pre_post_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, norm_struct, varargin{:});

size_dim2 = cellfun(@(x) size(x, 2), SW_xPlt.data);

%% Getting max. over bands.

load('seven_bands'), bands(3, 2) = 14; % bands(end, :) = [];
% bands = [1 4; 4 8; 8 14; 15 22; 23 30; 31 60; 60 100; 120 180];

[~, band_labels] = band_max(SW_xPlt.data{1}, frequencies, bands(1:end, :));

SW_Bands.data = cellfun(@(x) band_max(x, frequencies, bands(1:end, :)), SW_xPlt.data, 'UniformOutput', 0);

SW_Bands.meta.matrix_dim_1.name = '';
SW_Bands.meta.matrix_dim_1.values = band_labels;

load('M1_groups')

M1_increased_index{2} = M1_increased_index{2} & All_index{2};

M1_not_increased_index{2} = M1_not_increased_index{2} & All_index{2};

groups_plotted = {All_index, M1_increased_index, M1_not_increased_index}; % {All_index};

for group = 1:length(groups_plotted)

    group_name = [make_sliding_window_analysis_name([filename, groups_plotted{group}{1}, pd_label,...
    '_band', num2str(data_labels_struct.band_index), '_BP_', num2str(size(bands,1)), 'bands'], function_name,...
    window_time_cell, 2, varargin{:}), norm_label];

    SW_group = SW_Bands.packDim('Recording', 2);

    folders = SW_group.meta.matrix_dim_2.values;

    SW_group.data = cellfun(@(x) x(:, groups_plotted{group}{2}), SW_group.data, 'UniformOutput', 0);

    SW_group = SW_group.unpackDim(2);

    SW_group.axis(SW_group.findaxis('Recording')).values = ... % SW_xPlt = SW_xPlt.importAxisValues(SW_xPlt, 'Recording', {folders(groups_plotted{group}{2})});
        folders(groups_plotted{group}{2});

    % if ~any(size_dim2(:) > 1)

    if ~isempty(SW_group.findaxis('Window_Dim_1'))
        %% Plotting if there are sliding windows.

        close('all')

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

        function_handles = {@xp_tight_subplot_adaptive,@xp_compare_barplot_2D};
        function_arguments = {{},{@ttest, significance, [], 1}};
        dimensions = {{'Recording', 'Channel'},{'Period'}};
        recursivePlot(SW_WindowsPacked,function_handles,dimensions,function_arguments);

        save_all_figs([group_name, '_windows_by_subject_compared_p', num2str(significance)])

        close('all')

        % waitforbuttonpress

        %% % % Pre-post comparison of mean across windows, over recordings. % % %

        SW_RecordingsPacked = mean_over_axis(SW_group, 'Window_Dim_1');

        SW_RecordingsPacked = squeeze(SW_RecordingsPacked.packDim('Recording', 2));

        SW_RecordingsPacked.getaxisinfo

        dimensions = {{'Channel'},{'Period'}};
        recursivePlot(SW_RecordingsPacked,function_handles,dimensions,function_arguments);

        save_all_figs([group_name, '_windows_compared_p', num2str(significance)])

        save([group_name, '_recordingspacked.mat'], 'SW_RecordingsPacked')

    else
        %% Plotting if there are no sliding windows.

        close('all')

        %% % % Plots by recording. % % %

        SW_RecordingsPacked = squeeze(SW_group.packDim('Recording', 2));

        function_handles = {@xp_tight_subplot_adaptive,@xp_matrix_barplot};
        function_arguments = {{{'Period', 'Channel'}},{}};
        dimensions = {1:length(size(SW_RecordingsPacked)),0};
        recursivePlot(SW_RecordingsPacked,function_handles,dimensions,function_arguments);

        save_all_figs([group_name, '_recordings'])

        close('all')

        % waitforbuttonpress

        %% % % Pre-post comparisons over recordings. % % %

        function_handles = {@xp_tight_subplot_adaptive,@xp_compare_barplot_2D};
        function_arguments = {{},{@ttest, significance, [], 1}};
        dimensions = {{'Channel'},{'Period'}};
        recursivePlot(SW_RecordingsPacked,function_handles,dimensions,function_arguments);

        save_all_figs([group_name, '_recordings_compared_p', num2str(significance)])

        save([group_name, '_recordingspacked.mat'], 'SW_RecordingsPacked')

    end

    % else


end

end

function [bm, band_labels] = band_max(data, frequencies, bands)

no_bands = size(bands, 1);

if isnumeric(frequencies)

    bm = nan(no_bands, 1); band_labels = {};

    for b = 1:no_bands

        band_labels{b} = sprintf('%d-%d Hz', bands(b, :));

        band_indicator = frequencies >= bands(b, 1) & frequencies <= bands(b, 2);

        bm(b, 1) = max(data(band_indicator));

    end

elseif iscell(frequencies)

    freq_mat_1 = repmat(frequencies{1}', 1, length(frequencies{2}));

    freq_mat_2 = repmat(frequencies{2}, length(frequencies{1}), 1);

    bm = nan(no_bands); band_labels = {};

    bm_pair_index = 1; bm_indices = nan(no_bands, 2);

    for b2 = 1:no_bands

        band_indicator_2 = freq_mat_2 >= bands(b2, 1) & freq_mat_2 <= bands(b2, 2);

        for b1 = 1:no_bands

            band_indicator_1 = freq_mat_1 >= bands(b1, 1) & freq_mat_1 <= bands(b1, 2);

            band_indicator = band_indicator_1 & band_indicator_2;

            if ~isempty(data(band_indicator))

                bm(b1, b2) = max(data(band_indicator));

                band_labels{bm_pair_index} = sprintf('%d-%dx%d-%d', bands(b1, :), bands(b2, :));

                bm_pair_index = bm_pair_index + 1;

            end

        end

    end

    bm = bm(:);

    bm(isnan(bm)) = [];

end

end
