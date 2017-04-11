function SW_xPlt = PD_sliding_window_pre_post_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, norm, varargin)

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
  
        if isempty(SW_Baseline.findaxis('Period'))
            
            period_unpack_dim = SW_Baseline.lastNonSingletonDim + 1;
            
            SW_Baseline = SW_Baseline.unpackDim(period_unpack_dim, SW_xPlt.findaxis('Period'), 'Period', {'baseline'});
            
        end
        
        SW_Baseline = SW_Baseline.alignAxes(SW_xPlt);
        
        SW_BaselineMerged = SW_xPlt.merge(SW_Baseline);
        
        SW_BaselineMerged.getaxisinfo
        
        SW_xPlt = norm_axis_by_value(SW_BaselineMerged, 'Period', 'baseline');
        
        SW_xPlt = SW_xPlt.axissubset('Period', 'p');
        
        SW_xPlt.getaxisinfo
        
        % SW_xPlt.data = cellfun(@(x, y) 100*(x./y - 1), SW_xPlt.data, SW_Baseline.data, 'UniformOutput', 0);
        
    case 'shuffle'
    %% Normalization by shuffled surrogate data.

        SW_Shuffle = load([make_sliding_window_analysis_name([filename, '_pre_shuffles_post_shuffles_band',...
            num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        axes_info_struct = SW_Shuffle.axes_info_struct;
        
        SW_Shuffle = SW_Shuffle.SW_xPlt;
        
        if ~strcmp(function_name, 'PAC')
            
            SW_Shuffle = xp_abs(SW_Shuffle);
            
        end
        
        if ~isempty(SW_Shuffle.findaxis('Shuffles'))
            
            SW_ShuffleMean = mean_over_axis(SW_ShuffleMean, 'Shuffles');
            SW_ShuffleSTD = mean_over_axis(SW_ShuffleSTD, 'Shuffles', @nanstd, []);
            
        else
            
            [SW_ShuffleMean, SW_ShuffleSTD] = deal(SW_Shuffle);
            
            SW_ShuffleMean.data = cellfun(@(x) nanmean(x, axes_info_struct.odims + 1), SW_Shuffle.data, 'UniformOutput', false);
            SW_ShuffleSTD.data = cellfun(@(x) nanstd(x, [], axes_info_struct.odims + 1), SW_Shuffle.data, 'UniformOutput', false);
            
            SW_ShuffleMean.meta = rmfield(SW_ShuffleMean.meta, ['matrix_dim_', num2str(axes_info_struct.odims + 1)]);
            SW_ShuffleSTD.meta = rmfield(SW_ShuffleSTD.meta, ['matrix_dim_', num2str(axes_info_struct.odims + 1)]);
            
        end
        
        if ~isempty(SW_xPlt.findaxis('Window_Dim_1'))
            
            SW_ShuffleMean = SW_ShuffleMean.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
            SW_ShuffleSTD = SW_ShuffleSTD.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
            
        end
  
        shuffle_period_dim = SW_ShuffleMean.findaxis('Period');
        
        SW_ShuffleMean.axis(shuffle_period_dim).values = SW_xPlt.axis(SW_xPlt.findaxis('Period')).values;
        SW_ShuffleSTD.axis(shuffle_period_dim).values = SW_xPlt.axis(SW_xPlt.findaxis('Period')).values;
        
        data_type_axis = nDDictAxis;
        
        data_type_axis.name = 'Data_Type'; data_type_axis.values = {'Surrogate'};
        
        SW_ShuffleMean.axis(end + 1) = data_type_axis;
        SW_ShuffleSTD.axis(end + 1) = data_type_axis;
        
        data_type_axis = nDDictAxis;
        
        data_type_axis.name = 'Data_Type'; data_type_axis.values = {'Observation'};
        
        SW_xPlt.axis(end + 1) = data_type_axis;
        
        SW_ShuffleMean = SW_ShuffleMean.alignAxes(SW_xPlt);
        
        SW_ShuffleMeanMerged = SW_xPlt.merge(SW_ShuffleMean);
        
        SW_xPlt = norm_axis_by_value(SW_ShuffleMeanMerged, 'Data_Type', 'Surrogate', 'subtract');
        
        SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
        
        SW_ShuffleSTD = SW_ShuffleSTD.alignAxes(SW_xPlt);
        
        SW_ShuffleSTDMerged = SW_xPlt.merge(SW_ShuffleSTD);
        
        SW_xPlt = norm_axis_by_value(SW_ShuffleSTDMerged, 'Data_Type', 'Surrogate', 'divide');
        
        SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
        
        SW_xPlt.getaxisinfo
        
    case 'shuffle_subtracted'
    %% Normalization by shuffled surrogate data.

        SW_Shuffle = load([make_sliding_window_analysis_name([filename, '_pre_shuffles_post_shuffles_band',...
            num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        axes_info_struct = SW_Shuffle.axes_info_struct;
        
        SW_Shuffle = SW_Shuffle.SW_xPlt;
        
        if ~strcmp(function_name, 'PAC')
            
            SW_Shuffle = xp_abs(SW_Shuffle);
            
        end
        
        if ~isempty(SW_Shuffle.findaxis('Shuffles'))
            
            SW_ShuffleMean = mean_over_axis(SW_ShuffleMean, 'Shuffles');
            
        else
            
            SW_ShuffleMean = SW_Shuffle;
            
            SW_ShuffleMean.data = cellfun(@(x) nanmean(x, axes_info_struct.odims + 1), SW_Shuffle.data, 'UniformOutput', false);
            
            SW_ShuffleMean.meta = rmfield(SW_ShuffleMean.meta, ['matrix_dim_', num2str(axes_info_struct.odims + 1)]);
            
        end
        
        if ~isempty(SW_xPlt.findaxis('Window_Dim_1'))
            
            % wd1_axis = SW_xPlt.findaxis('Window_Dim_1');
            % window_pack_dim = SW_xPlt.lastNonSingletonDim + 1;
            % SW_xPlt = SW_xPlt.packDim('Window_Dim_1', SW_xPlt.lastNonSingletonDim + 1);
            % SW_xPlt = SW_xPlt.unpackDim(window_pack_dim, wd1_axis);
            
            SW_ShuffleMean = SW_ShuffleMean.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
            
        end
  
        shuffle_period_dim = SW_ShuffleMean.findaxis('Period');
        data_period_dim = SW_xPlt.findaxis('Period');
        
        if ~isempty(data_period_dim) && ~isempty(shuffle_period_dim)
            
            SW_ShuffleMean.axis(shuffle_period_dim).values = SW_xPlt.axis(SW_xPlt.findaxis('Period')).values;
            
        else
            
            if isempty(shuffle_period_dim)
                
                period_unpack_dim = SW_ShuffleMean.lastNonSingletonDim + 1;
                
                SW_ShuffleMean = SW_ShuffleMean.unpackDim(period_unpack_dim, data_period_dim, 'Period', pd_names);
                
            end
            
            if isempty(data_period_dim)
                
                period_unpack_dim = SW_xPlt.lastNonSingletonDim + 1;
                
                SW_xPlt = SW_xPlt.unpackDim(period_unpack_dim, shuffle_period_dim, 'Period', pd_names);
                
            end
            
        end
        
        data_type_axis = nDDictAxis;
        
        data_type_axis.name = 'Data_Type'; data_type_axis.values = {'Surrogate'};
        
        SW_ShuffleMean.axis(end + 1) = data_type_axis;
        
        data_type_axis = nDDictAxis;
        
        data_type_axis.name = 'Data_Type'; data_type_axis.values = {'Observation'};
        
        SW_xPlt.axis(end + 1) = data_type_axis;
        
        SW_ShuffleMean = SW_ShuffleMean.alignAxes(SW_xPlt);
        
        SW_ShuffleMeanMerged = SW_xPlt.merge(SW_ShuffleMean);
        
        SW_xPlt = norm_axis_by_value(SW_ShuffleMeanMerged, 'Data_Type', 'Surrogate', 'subtract');
        
        SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
        
        SW_xPlt.getaxisinfo
        
    case 'shuffle_baseline'
    %% Normalization by shuffled data followed by normalization by (shuffle-normalized) baseline data.

        SW_xPlt = PD_sliding_window_pre_post_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, 'shuffle', varargin{:});
        
        SW_xPlt.getaxisinfo

        SW_Baseline = load([make_sliding_window_analysis_name([filename, '_baseline_band',...
            num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        SW_Baseline = SW_Baseline.SW_xPlt;
        
        SW_Baseline = PD_sliding_window_normalize(SW_Baseline, function_name, sliding_window_cell, data_labels_struct, filename, {'baseline'}, 'shuffle', {}, varargin{:});
        
        if ~strcmp(function_name, 'PAC')
            
            SW_Baseline = xp_abs(SW_Baseline);
        
        end
        
        if ~isempty(SW_Baseline.findaxis('Window_Dim_1'))
            
            SW_Baseline = mean_over_axis(SW_Baseline, 'Window_Dim_1');
            
        end
        
        if ~isempty(SW_xPlt.findaxis('Window_Dim_1'))
            
            SW_Baseline = SW_Baseline.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
            
        end
        
        SW_Baseline = SW_Baseline.alignAxes(SW_xPlt);
        
        SW_BaselineMerged = SW_xPlt.merge(SW_Baseline);
        
        % SW_BaselineMerged.getaxisinfo
        
        SW_xPlt = norm_axis_by_value(SW_BaselineMerged, 'Period', 'baseline', 'subtract');
        
        SW_xPlt = SW_xPlt.axissubset('Period', 'p');
        
        SW_xPlt.getaxisinfo
        
    case 'baseline_shuffle'
    %% Normalization by baseline data followed by normalization by (baseline-normalized) shuffle data.

        SW_xPlt = PD_sliding_window_pre_post_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, 'baseline', varargin{:});
        
        SW_xPlt.getaxisinfo

        SW_Shuffle = load([make_sliding_window_analysis_name([filename, '_pre_shuffles_post_shuffles_band',...
            num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        axes_info_struct = SW_Shuffle.axes_info_struct;
        
        SW_Shuffle = SW_Shuffle.SW_xPlt;
        
        SW_Shuffle = SW_Shuffle.unpackDim(axes_info_struct.odims + 1, [], 'Window_Dim_1');
        
        SW_Shuffle = PD_sliding_window_normalize(SW_Shuffle, function_name, sliding_window_cell, data_labels_struct, filename, {'pre_shuffles', 'post_shuffles'}, 'baseline', {'baseline_shuffles'}, varargin{:});
        
        if ~strcmp(function_name, 'PAC')
            
            SW_Shuffle = xp_abs(SW_Shuffle);
            
        end
        
        SW_ShuffleMean = mean_over_axis(SW_Shuffle, 'Window_Dim_1');
        SW_ShuffleSTD = mean_over_axis(SW_Shuffle, 'Window_Dim_1', @nanstd, []);
        
        if ~isempty(SW_xPlt.findaxis('Window_Dim_1'))
            
            SW_ShuffleMean = SW_ShuffleMean.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
            SW_ShuffleSTD = SW_ShuffleSTD.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
            
        end
  
        shuffle_period_dim = SW_ShuffleMean.findaxis('Period');
        
        SW_ShuffleMean.axis(shuffle_period_dim).values = SW_xPlt.axis(SW_xPlt.findaxis('Period')).values;
        SW_ShuffleSTD.axis(shuffle_period_dim).values = SW_xPlt.axis(SW_xPlt.findaxis('Period')).values;
        
        data_type_axis = nDDictAxis;
        
        data_type_axis.name = 'Data_Type'; data_type_axis.values = {'Surrogate'};
        
        SW_ShuffleMean.axis(end + 1) = data_type_axis;
        SW_ShuffleSTD.axis(end + 1) = data_type_axis;
        
        data_type_axis = nDDictAxis;
        
        data_type_axis.name = 'Data_Type'; data_type_axis.values = {'Observation'};
        
        SW_xPlt.axis(end + 1) = data_type_axis;
        
        SW_ShuffleMean = SW_ShuffleMean.alignAxes(SW_xPlt);
        
        SW_ShuffleMeanMerged = SW_xPlt.merge(SW_ShuffleMean);
        
        SW_xPlt = norm_axis_by_value(SW_ShuffleMeanMerged, 'Data_Type', 'Surrogate', 'subtract');
        
        SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
        
        SW_ShuffleSTD = SW_ShuffleSTD.alignAxes(SW_xPlt);
        
        SW_ShuffleSTDMerged = SW_xPlt.merge(SW_ShuffleSTD);
        
        SW_xPlt = norm_axis_by_value(SW_ShuffleSTDMerged, 'Data_Type', 'Surrogate', 'divide');
        
        SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
        
        SW_xPlt.getaxisinfo
            
end