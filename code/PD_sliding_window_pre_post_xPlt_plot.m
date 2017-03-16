function PD_sliding_window_pre_post_xPlt_plot(function_name, sliding_window_cell, data_labels_struct, filename, varargin)
    
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

pd_names = {'Pre-Infusion', 'Post-Infusion'};

load([make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index)], function_name,...
    window_time_cell, 2, varargin{:}), '_xPlt.mat'])

SW_xPlt = abs(SW_xPlt);

SW_xPlt.getaxisinfo

if ~isempty(SW_xPlt.axis.findAxes('Window_Dim_1'))

    SW_Windows = squeeze(SW_xPlt.packDim('Window_Dim_1', 2));
    
    SW_Windows.getaxisinfo
    
    function_handles = {@xp_subplot_grid_adaptive,@xp_matrix};
    function_arguments = {{{'Recording', 'Period', 'To', 'From'}},{}};
    dimensions = {1:length(size(SW_Windows)),0};
    recursivePlot(SW_Windows,function_handles,dimensions,function_arguments);
    
end

SW_Recordings = squeeze(SW_xPlt.packDim('Recording', 2));

function_handles = {@xp_subplot_grid_adaptive,@xp_matrix_basicplot};
function_arguments = {{{'Windows', 'Period', 'To', 'From'}},{}};
dimensions = {1:length(size(SW_Recordings)),0};
recursivePlot(SW_Recordings,function_handles,dimensions,function_arguments);

end

