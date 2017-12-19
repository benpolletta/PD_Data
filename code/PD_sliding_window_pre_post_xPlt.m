function SW_xPlt = PD_sliding_window_pre_post_xPlt(function_name, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, varargin)
    
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
    '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.bands), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}))

SW_size = size(SW);

SW_xPlt = xPlt;

SW_xPlt = SW_xPlt.importData({SW});

no_windows = cellfun(@(x) length(x), window_time);

no_windows(no_windows == 1) = [];

wdims = length(no_windows);
    
output_size(output_size == 1) = [];

odims = length(output_size);

axes_info_struct = get_axes_info(function_name,...
    sliding_window_cell, data_labels_struct, no_windows, output_size, varargin{:});

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

for wdim = 1:wdims
    
    SW_xPlt = unpackDim(SW_xPlt, length(SW_size) - dims_from_last, 1,...
        axes_info_struct.window_names{wdims - wdim + 1}, axes_info_struct.window_values{wdims - wdim + 1});
    
    dims_from_last = dims_from_last + 1;
    
end
    
%% Unpacking output dimensions.

for odim = 1:(odims - 1)
    
    SW_xPlt = unpackDim(SW_xPlt, length(SW_size) - dims_from_last, 1,...
        axes_info_struct.output_names{odims - odim + 1}, axes_info_struct.output_values{odims - odim + 1});
    
    dims_from_last = dims_from_last + 1;
    
end

remaining_axis = nDDictAxis;

remaining_axis.name = axes_info_struct.output_names{1};
remaining_axis.values = axes_info_struct.output_values{1};

meta = SW_xPlt.meta;

meta.matrix_dim_1 = remaining_axis;

SW_xPlt = importMeta(SW_xPlt, meta);

SW_xPlt = squeeze(SW_xPlt);

save([make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.bands), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}), '_xPlt.mat'], 'SW_xPlt')

end

function axes_info_struct = get_axes_info(function_name,...
    sliding_window_cell, data_labels_struct, window_size, output_size, varargin)

axes_info_struct.output_names{1} = 'Freq. (Hz)';

axes_info_struct.output_values{1} = data_labels_struct.sampling_freq{1}*...
    ((1:output_size(1)) - 1)/(2*(output_size(1) - 1));

for odim = 2:length(output_size)
    
    axes_info_struct.output_names{odim} = ['Output_Dim_' num2str(odim)];
    
    axes_info_struct.output_values{odim} = 1:output_size(odim);
    
end

for wdim = 1:length(window_size)
    
    axes_info_struct.window_names{wdim} = ['Window_Dim_' num2str(wdim)];
    
    axes_info_struct.window_values{wdim} = 1:window_size(wdim);
    
end

switch function_name
    
    case 'mvgc_analysis'
        
        switch varargin{end}
            
            case 1
                
                axes_info_struct.output_names([2 3]) = {'Channel To', 'Channel From'};
                
                axes_info_struct.output_values([2 3]) = deal({{'Striatum', 'Motor Ctx.'}});
                
            case 3
                
                if sliding_window_cell{2}(1) == 1
                    
                    if sliding_window_cell{1}(1) == 150*500
                        
                        axes_info_struct.window_names{1} = 'Channel';
                        
                        axes_info_struct.window_values{1} = {'Striatum', 'Motor Ctx.'};
                        
                    else
                        
                        axes_info_struct.window_names{2} = 'Channel';
                        
                        axes_info_struct.window_values{2} = {'Striatum', 'Motor Ctx.'};
                        
                    end
                    
                elseif sliding_window_cell{2}(1) == 2
                    
                    axes_info_struct.output_names([2 3]) = {'Channel 1', 'Channel 2'};
                
                    axes_info_struct.output_values([2 3]) = deal({{'Striatum', 'Motor Ctx.'}});
                    
                end
                
        end
        
    case 'pmtm'
        
        if sliding_window_cell{1}(1) == 150*500
            
            axes_info_struct.window_names{1} = 'Channel';
            
            axes_info_struct.window_values{1} = {'Striatum', 'Motor Ctx.'};
            
        else
            
            axes_info_struct.window_names{2} = 'Channel';
            
            axes_info_struct.window_values{2} = {'Striatum', 'Motor Ctx.'};
            
        end
        
    case 'pmtm_detrend'
        
        if sliding_window_cell{1}(1) == 150*500
            
            axes_info_struct.window_names{1} = 'Channel';
            
            axes_info_struct.window_values{1} = {'Striatum', 'Motor Ctx.'};
            
        else
            
            axes_info_struct.window_names{2} = 'Channel';
            
            axes_info_struct.window_values{2} = {'Striatum', 'Motor Ctx.'};
            
        end
        
end

end

