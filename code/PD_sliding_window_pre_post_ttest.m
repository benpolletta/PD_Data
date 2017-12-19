function p_vals = PD_sliding_window_pre_post_ttest(filename, comparison_struct, function_name, sliding_window_cell, data_labels_struct, varargin)

% Performs t-test on carbachol data at times of highest
% striatal beta band density, pre- and post-infusion.
% INPUTS:
%     filename (string): prefix for analysis, e.g. 'STR_w_M1'.
%     comparison_struct: structure which can be initialized with
%       init_comparisons.m, containing fields:
%       comparison_struct.comparison_indices (N_comparisons x 2 cell of N_dims x 1 cells): cell
%           containing indices of sliding window data to be compared against
%           each other.
%       comparison_size (N_dims x 1 vector): how comparisons should be
%           resized (i.e., what dimensions are looped over in
%           "linearized" cell of comparisons comparison_struct.comparison_indices).
%       comparison_name (string): name to save the comparison.
%     function_name (string or function handle): (name of) function
%       performed on each window.
%     sliding_window_cell (2 x 1 cell of 2 x 1 arrays): cell
%       containing sliding window length and step length (in indices) for the
%       two dimensions of the carbachol data.
%     data_labels_struct (structure): can be initialized by
%       init_data_labels.m, contains fields: BP_suffix, peak_suffix,
%       data_suffix, epoch_secs, pd_handle, data_field, band_index,
%       sampling_freq.

function_name = get_fname(function_name);

window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

analysis_name = make_sliding_window_analysis_name([filename, '_pre_post',...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.bands), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:});

SW = load(analysis_name); SW = SW.SW;

number_comparisons = length(comparison_struct.comparison_indices);

p_vals = nan(number_comparisons, 2);
    
parfor comparison = 1:number_comparisons
    
    temp = nan(1, 2);
    
    [~, temp(1)] = ttest(SW(comparison_struct.comparison_indices{comparison, 1}{:}),...
        SW(comparison_struct.comparison_indices{comparison, 2}{:}), 'tail', 'left');
    
    [~, temp(2)] = ttest(SW(comparison_struct.comparison_indices{comparison, 1}{:}),...
        SW(comparison_struct.comparison_indices{comparison, 2}{:}), 'tail', 'right');
    
    p_vals(comparison, :) = temp;
    
end

p_vals = reshape(p_vals, [comparison_struct.comparison_size, 2]);

mat_name = make_sliding_window_analysis_name([filename, comparison_struct.comparison_name,...
    '_pre_post_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.bands), 7, 'back')],...
    function_name, window_time_cell, 2, varargin{:});

save([mat_name, '_ttest.mat'], 'p_vals')

