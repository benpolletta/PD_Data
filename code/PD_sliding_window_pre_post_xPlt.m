function SW_xPlt = PD_sliding_window_pre_post_xPlt(function_name, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, output_struct, varargin)
    
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

load(make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index)], function_name,...
    window_time_cell, 2, varargin{:}))

SW_size = size(SW);

SW_xPlt = xPlt;

SW_xPlt = SW_xPlt.importData({SW});

% SW_xPlt = SW_xPlt.fixAxes;

dims_from_last = 0;

%% Unpacking period (last dimension).

pd_names = {'Pre-Infusion', 'Post-Infusion'};

SW_xPlt = unpackDim(SW_xPlt, length(SW_size) - dims_from_last, 1, 'Period', pd_names);

dims_from_last = dims_from_last + 1;

%% Unpacking recordings (second-to-last dimension).

folder_index = 0;

for s = 1:length(subjects_mat_cell)
    
    load(subjects_mat_cell{s})
    
    All_folders(folder_index + (1:length(folders))) = folders;
    
    folder_index = folder_index + length(folders);
    
end

SW_xPlt = unpackDim(SW_xPlt, length(SW_size) - dims_from_last, 1, 'Recording', All_folders);

dims_from_last = dims_from_last + 1;

%% Unpacking window dimensions.

total_windows = cellfun(@(x) length(x), window_time);

total_windows(total_windows == 1) = [];

wdims = length(total_windows);

for wdim = 1:wdims
    
    SW_xPlt = unpackDim(SW_xPlt, length(SW_size) - dims_from_last, 1, ['Window_Dim_' num2str(wdims - wdim + 1)], window_time{wdims - wdim + 1});
    
    dims_from_last = dims_from_last + 1;
    
end
    
%% Unpacking output dimensions.

if output_struct.unpack_flag
    
    output_size(output_size == 1) = [];
    
    output_struct = init_output_axis(output_size, output_struct);
    
    odims = length(output_size);
    
    for odim = 1:(odims - 1)
        
        SW_xPlt = unpackDim(SW_xPlt, length(SW_size) - dims_from_last, 1, output_struct.output_names{odims - odim + 1}, output_struct.output_values{odims - odim + 1});
        
        dims_from_last = dims_from_last + 1;
        
    end
    
    remaining_axis = nDDictAxis;
    
    remaining_axis.name = output_struct.output_names{1};
    remaining_axis.values = output_struct.output_values{1};
    
    meta = SW_xPlt.meta;
    
    meta.matrix_dim_1 = remaining_axis;
    
    SW_xPlt = importMeta(SW_xPlt, meta);
    
end

SW_xPlt = squeeze(SW_xPlt);

save([make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index)], function_name,...
    window_time_cell, 2, varargin{:}), '_xPlt.mat'], 'SW_xPlt')

end

function output_struct = init_output_axis(output_size, output_struct)
    
    if ~isfield(output_struct, 'output_names'), output_struct.output_names = []; end
    
    if isempty(output_struct.output_names)
       
        for odim = 1:length(output_size)
            
           output_names{odim} = ['Output_Dim_' num2str(odim)];
            
        end
        
        output_struct.output_names = output_names;
        
    end

    if ~isfield(output_struct, 'output_values'), output_struct.output_values = []; end
    
    if isempty(output_struct.output_values)
       
        for odim = 1:length(output_size)
            
           output_values{odim} = 1:output_size(odim);
            
        end
        
        output_struct.output_values = output_values;
        
    end

end

