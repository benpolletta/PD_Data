function SW_xPlt = PD_sliding_window_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, pd_names, norm, varargin)

function_name = get_fname(function_name);
    
window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

switch norm
    
    case ''   
    %% No normalization.
    
    case 'totalpower'
    %% Normalization by total power.
        
        SW_xPlt.data = cellfun(@(x) x/sum(x), SW_xPlt.data, 'UniformOutput', 0);
    
    case 'frequency'
    %% Normalization by frequency (to remove 1/f).
        
        freq_mult = diag(frequencies.^(2/3));
        
        SW_xPlt.data = cellfun(@(x) freq_mult*x, SW_xPlt.data, 'UniformOutput', 0);
        
    case 'baseline'
    %% Normalization by mean baseline value.

        SW_Baseline = load([make_sliding_window_analysis_name([filename, '_baseline_band',...
            num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        SW_Baseline = SW_Baseline.SW_xPlt;
        
        if ~strcmp(function_name, 'PAC')
            
            SW_Baseline = xp_abs(SW_Baseline);
        
        end
        
        if ~isempty(SW_Baseline.findaxis('Window_Dim_1'))
            
            SW_Baseline = mean_over_axis(SW_Baseline, 'Window_Dim_1');
            
        end
        
        % SW_Baseline = SW_Baseline.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Period')).values, 'Period');
        
        if ~isempty(SW_xPlt.findaxis('Window_Dim_1'))
            
            SW_Baseline = SW_Baseline.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
            
        end
  
        period_unpack_dim = SW_Baseline.lastNonSingletonDim + 1;
        
        SW_Baseline = SW_Baseline.unpackDim(period_unpack_dim, SW_xPlt.findaxis('Period'), 'Period', {'baseline'});
        
        SW_Baseline.getaxisinfo
        
        SW_BaselineMerged = SW_xPlt.merge(SW_Baseline);
        
        SW_BaselineMerged.getaxisinfo
        
        SW_xPlt = norm_axis_by_value(SW_BaselineMerged, 'Period', 'baseline');
        
        SW_xPlt = SW_xPlt.axissubset('Period', 'p');
        
        SW_xPlt.getaxisinfo
        
        % SW_xPlt.data = cellfun(@(x, y) 100*(x./y - 1), SW_xPlt.data, SW_Baseline.data, 'UniformOutput', 0);
        
    case 'shuffle'
    %% Normalization by shuffled surrogate data.
    
        shuffle_label = '';
    
        for pd = 1:length(pd_names)
        
            shuffle_label = [shuffle_label, '_', pd_names{pd}, '_shuffles'];
            
        end

        SW_Shuffle = load([make_sliding_window_analysis_name([filename, shuffle_label,...
            num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        axes_info_struct = SW_Shuffle.axes_info_struct;
        
        SW_Shuffle = SW_Shuffle.SW_xPlt;
        
        if ~strcmp(function_name, 'PAC')
            
            SW_Shuffle = xp_abs(SW_Shuffle);
            
        end
        
        SW_Shuffle.data = cellfun(@(x) nanmean(x, axes_info_struct.leave_packed_odim + 1), SW_Shuffle.data, 'UniformOutput', false);
        
        SW_Shuffle.meta = rmfield(SW_Shuffle.meta, ['matrix_dim_', num2str(axes_info_struct.leave_packed_odim + 1)]);
        
        if ~isempty(SW_Shuffle.findaxis('Window_Dim_1'))
            
            SW_Shuffle = mean_over_axis(SW_Shuffle, 'Window_Dim_1');
            
        end
        
        if ~isempty(SW_xPlt.findaxis('Window_Dim_1'))
            
            SW_Shuffle = SW_Shuffle.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
            
        end
  
        shuffle_period_dim = SW_Shuffle.findaxis('Period');
        
        SW_Shuffle.axis(shuffle_period_dim).values = SW_xPlt.axis(SW_xPlt.findaxis('Period')).values;
        
        data_type_axis = nDDictAxis;
        
        data_type_axis.name = 'Data_Type'; data_type_axis.values = {'Surrogate'};
        
        SW_Shuffle.axis(end + 1) = data_type_axis;
        
        data_type_axis = nDDictAxis;
        
        data_type_axis.name = 'Data_Type'; data_type_axis.values = {'Observation'};
        
        SW_xPlt.axis(end + 1) = data_type_axis;
        
        SW_xPlt.getaxisinfo
        
        SW_Shuffle.getaxisinfo
        
        SW_ShuffleMerged = SW_xPlt.merge(SW_Shuffle);
        
        SW_ShuffleMerged.getaxisinfo
        
        SW_xPlt = norm_axis_by_value(SW_ShuffleMerged, 'Data_Type', 'Surrogate');
        
        SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
        
        SW_xPlt.getaxisinfo
        
    case 'shuffle_baseline'
    %% Normalization by shuffled data followed by normalization by (shuffle-normalized) baseline data.

        SW_xPlt = PD_sliding_window_pre_post_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, {'pre', 'post'}, 'shuffle', varargin{:});
        
        SW_xPlt.getaxisinfo

        SW_Baseline = load([make_sliding_window_analysis_name([filename, '_baseline_band',...
            num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        SW_Baseline = SW_Baseline.SW_xPlt;
        
        SW_Baseline = PD_sliding_window_pre_post_normalize(SW_Baseline, function_name, sliding_window_cell, data_labels_struct, filename, {'baseline'}, 'shuffle', varargin{:});
        
        if ~strcmp(function_name, 'PAC')
            
            SW_Baseline = xp_abs(SW_Baseline);
        
        end
        
        if ~isempty(SW_Baseline.findaxis('Window_Dim_1'))
            
            SW_Baseline = mean_over_axis(SW_Baseline, 'Window_Dim_1');
            
        end
        
        if ~isempty(SW_xPlt.findaxis('Window_Dim_1'))
            
            SW_Baseline = SW_Baseline.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
            
        end
  
        period_unpack_dim = SW_Baseline.lastNonSingletonDim + 1;
        
        SW_Baseline = SW_Baseline.unpackDim(period_unpack_dim, SW_xPlt.findaxis('Period'), 'Period', {'baseline'});
        
        SW_Baseline.getaxisinfo
        
        SW_BaselineMerged = SW_xPlt.merge(SW_Baseline);
        
        SW_BaselineMerged.getaxisinfo
        
        SW_xPlt = norm_axis_by_value(SW_BaselineMerged, 'Period', 'baseline');
        
        SW_xPlt = SW_xPlt.axissubset('Period', 'p');
        
        SW_xPlt.getaxisinfo
            
end