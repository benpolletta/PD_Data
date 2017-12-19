function [sw_xPlt, axes_info_struct] = PD_sliding_window_xPlt(function_name, sliding_window_cell, subjects_mat_struct, data_labels_struct, filename, pd_names, varargin)
    
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

no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods
    
    pd_label = [pd_label, '_', pd_names{period}];
    
end

load([make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.bands), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}), '.mat'])

sw_size = size(sw);

sw_xPlt = xPlt;

sw_xPlt = sw_xPlt.importData({sw});

meta = sw_xPlt.meta;

if ~exist('no_windows', 'var'), no_windows = cellfun(@(x) length(x), window_time); end

no_windows(no_windows == 1) = [];

wdims = length(no_windows);
    
output_size(output_size == 1) = [];

odims = length(output_size);

axes_info_struct = get_axes_info(function_name, sliding_window_cell,...
    subjects_mat_struct.chan_labels, data_labels_struct, pd_names, no_windows, output_size, varargin{:});

dims_from_last = 0;

%% Unpacking window dimensions.

leave_packed_wdim = [];

for wdim = wdims:-1:1
    
    if no_windows(wdim) < 100
        
        sw_xPlt = unpackDim(sw_xPlt, length(sw_size) - dims_from_last, 1,...
            axes_info_struct.window_names{wdim}, axes_info_struct.window_values{wdim});
        
    else
        
        remaining_axis = nDDictAxis;
        
        remaining_axis.name = axes_info_struct.window_names{wdim};
        remaining_axis.values = axes_info_struct.window_values{wdim};
        
        meta.(['matrix_dim_', num2str(odims + wdim)]) = remaining_axis;
        
        leave_packed_wdim(end + 1) = wdim;
        
    end
    
    dims_from_last = dims_from_last + 1;
    
end
    
if ~isempty(leave_packed_wdim), axes_info_struct.leave_packed_wdim = leave_packed_wdim; end

%% Unpacking output dimensions.

for odim = 1:(odims - axes_info_struct.leave_packed_odim)
    
    sw_xPlt = unpackDim(sw_xPlt, length(sw_size) - dims_from_last, 1,...
        axes_info_struct.output_names{odims - odim + 1}, axes_info_struct.output_values{odims - odim + 1});
    
    dims_from_last = dims_from_last + 1;
    
end

for odim = 1:axes_info_struct.leave_packed_odim
    
    remaining_axis = nDDictAxis;
    
    remaining_axis.name = axes_info_struct.output_names{odim};
    remaining_axis.values = axes_info_struct.output_values{odim};
    
    meta.(['matrix_dim_', num2str(odim)]) = remaining_axis;
    
end

sw_xPlt = importMeta(sw_xPlt, meta);

sw_xPlt = squeeze(sw_xPlt);

save([make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index), make_label('bands', length(data_labels_struct.bands), 7, 'back')], function_name,...
    window_time_cell, 2, varargin{:}), '_xPlt.mat'], 'sw_xPlt', 'axes_info_struct')

end

function axes_info_struct = get_axes_info(function_name, sliding_window_cell,...
    chan_labels, data_labels_struct, pd_names, window_size, output_size, varargin)

axes_info_struct.odims = length(output_size); axes_info_struct.wdims = length(window_size);

for odim = 1:length(output_size)
    
    axes_info_struct.output_names{odim} = ['Output_Dim_' num2str(odim)];
    
    axes_info_struct.output_values{odim} = 1:output_size(odim);
    
end

if strcmp(function_name, 'PAC')
    
    axes_info_struct.output_names{1} = 'Amp. Freq. (Hz)';
    
    axes_info_struct.output_values{1} = varargin{3};
    
    axes_info_struct.output_names{2} = 'Phase Freq. (Hz)';
    
    axes_info_struct.output_values{2} = varargin{2};
    
    axes_info_struct.leave_packed_odim = 2;
    
else
    
    axes_info_struct.output_names{1} = 'Freq. (Hz)';
    
    axes_info_struct.output_values{1} = data_labels_struct.sampling_freq{1}*...
        (0:(output_size(1) - 1))/(2*(output_size(1) - 1));
    
    axes_info_struct.leave_packed_odim = 1;
    
end

for wdim = 1:length(window_size)
    
    axes_info_struct.window_names{wdim} = ['Window_Dim_' num2str(wdim)];
    
    axes_info_struct.window_values{wdim} = 1:window_size(wdim);
    
end

if ~isempty(regexp(pd_names{:}, 'shuffle', 'once'))

    axes_info_struct.window_names{1} = 'Shuffles';
    
end

switch function_name
    
    case 'mvgc_analysis'
        
        switch varargin{end}
            
            case 1
                
                axes_info_struct.output_names([2 3]) = {'Channel To', 'Channel From'};
                
                axes_info_struct.output_values([2 3]) = deal({chan_labels});
                
            case 3
                
                if sliding_window_cell{2}(1) == 1
                    
                    if sliding_window_cell{1}(1) == 150*500 && ~strcmp(pd_names, 'baseline')
                        
                        axes_info_struct.window_names{1} = 'Channel';
                        
                        axes_info_struct.window_values{1} = chan_labels;
                        
                    else
                        
                        axes_info_struct.window_names{2} = 'Channel';
                        
                        axes_info_struct.window_values{2} = chan_labels;
                        
                    end
                    
                elseif sliding_window_cell{2}(1) == 2
                    
                    axes_info_struct.output_names([2 3]) = {'Channel 1', 'Channel 2'};
                
                    axes_info_struct.output_values([2 3]) = deal({chan_labels});
                    
                end
                
        end
        
    case 'pmtm'
        
        if sliding_window_cell{1}(1) == 150*500 && ~strcmp(pd_names, 'baseline')
            
            axes_info_struct.window_names{1} = 'Channel';
            
            axes_info_struct.window_values{1} = chan_labels;
            
        else
            
            axes_info_struct.window_names{2} = 'Channel';
            
            axes_info_struct.window_values{2} = chan_labels;
            
        end
        
    case 'pmtm_detrend'
        
        if sliding_window_cell{1}(1) == 150*500 && ~strcmp(pd_names, 'baseline')
            
            axes_info_struct.window_names{1} = 'Channel';
            
            axes_info_struct.window_values{1} = chan_labels;
            
        else
            
            axes_info_struct.window_names{2} = 'Channel';
            
            axes_info_struct.window_values{2} = chan_labels;
            
        end
        
    case 'PAC'
        
        if sliding_window_cell{1}(1) == 150*500 && ~strcmp(pd_names, 'baseline')
            
            axes_info_struct.window_names{1} = 'Channel';
            
            axes_info_struct.window_values{1} = chan_labels;
            
        else
            
            axes_info_struct.window_names{2} = 'Channel';
            
            axes_info_struct.window_values{2} = chan_labels;
            
        end
        
end

end

