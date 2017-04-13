function SW_xPlt = PD_sliding_window_pre_post_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, norm_struct, varargin)

function_name = get_fname(function_name);
    
window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

switch norm_struct.who
    
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

        SW_xPlt = PD_sliding_window_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, {'pre', 'post'}, norm_struct, {'baseline'}, varargin{:});

        SW_xPlt.getaxisinfo
        
        % SW_xPlt.data = cellfun(@(x, y) 100*(x./y - 1), SW_xPlt.data, SW_Baseline.data, 'UniformOutput', 0);
        
    case 'shuffle'
    %% Normalization by shuffled surrogate data.

        SW_xPlt = PD_sliding_window_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, {'pre', 'post'}, norm_struct, {'pre_shuffles', 'post_shuffles'}, varargin{:});
        
        SW_xPlt.getaxisinfo
        
    case 'shuffle_baseline'
    %% Normalization by shuffled data followed by normalization by (shuffle-normalized) baseline data.
        
        shuffle_prepost_struct = struct('who', 'shuffle', 'how', norm_struct.how{1});

        SW_xPlt = PD_sliding_window_pre_post_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, shuffle_prepost_struct, varargin{:});
        
        SW_xPlt.getaxisinfo

        SW_Baseline_struct = load([make_sliding_window_analysis_name([filename, '_baseline_band',...
            num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        shuffle_baseline_struct = struct('who', 'shuffle', 'how', norm_struct.how{1});
        
        SW_Baseline_struct.SW_xPlt = PD_sliding_window_normalize(SW_Baseline_struct.SW_xPlt,...
            function_name, sliding_window_cell, data_labels_struct, filename, {'baseline'}, shuffle_baseline_struct, {}, varargin{:});
        
        collapse_struct = struct('function', @nanmean, 'varargin', {});
        
        SW_Baseline = PD_xPlt_prep_for_norm(SW_Baseline_struct, SW_xPlt, function_name, {'baseline'}, collapse_struct);
        
        SW_BaselineMerged = SW_xPlt.merge(SW_Baseline);
        
        switch norm_struct.how{2}
            
            case ''
        
                SW_xPlt = norm_axis_by_value(SW_BaselineMerged, 'Period', 'baseline');
                
            case 'subtract'
        
                SW_xPlt = norm_axis_by_value(SW_BaselineMerged, 'Period', 'baseline', 'subtract');
                
            case 'zscore'
        
                SW_xPlt = norm_axis_by_value(SW_BaselineMerged, 'Period', 'baseline', 'subtract');
        
                SW_xPlt = SW_xPlt.axissubset('Period', 'p');
                
                collapse_struct = struct('function', @nanstd, 'varargin', []);
                
                SW_Baseline = PD_xPlt_prep_for_norm(SW_Baseline_struct, SW_xPlt, function_name, norm_pd_names, collapse_struct);
                
                SW_BaselineMerged = SW_xPlt.merge(SW_Baseline);
                
                SW_xPlt = norm_axis_by_value(SW_BaselineMerged, 'Period', 'baseline', 'divide');
                
        end
        
        SW_xPlt = SW_xPlt.axissubset('Period', 'p');
        
        SW_xPlt.getaxisinfo
        
    case 'baseline_shuffle'
    %% Normalization by baseline data followed by normalization by (baseline-normalized) shuffle data.
        
        baseline_prepost_struct = struct('who', 'baseline', 'how', norm_struct.how{1});

        SW_xPlt = PD_sliding_window_pre_post_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, baseline_prepost_struct, varargin{:});
        
        SW_xPlt.getaxisinfo

        SW_Shuffle_struct = load([make_sliding_window_analysis_name([filename, '_pre_shuffles_post_shuffles_band',...
            num2str(data_labels_struct.band_index)], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        baseline_shuffle_struct = struct('who', 'shuffle', 'how', norm_struct.how{1});
        
        SW_Shuffle_struct.SW_xPlt = PD_sliding_window_normalize(SW_Shuffle_struct.SW_xPlt,...
            function_name, sliding_window_cell, data_labels_struct, filename, {'pre_shuffles', 'post_shuffles'}, baseline_shuffle_struct, {'baseline_shuffles'}, varargin{:});
        
        collapse_struct = struct('function', @nanmean, 'varargin', {});
        
        SW_Shuffle = PD_xPlt_prep_for_norm(SW_Shuffle_struct, SW_xPlt, function_name, {'baseline'}, collapse_struct);
        
        surrogate_data_type_axis = nDDictAxis;
        
        surrogate_data_type_axis.name = 'Data_Type'; surrogate_data_type_axis.values = {'Surrogate'};
        
        SW_Shuffle.axis(end + 1) = surrogate_data_type_axis;
        
        observed_data_type_axis = nDDictAxis;
        
        observed_data_type_axis.name = 'Data_Type'; observed_data_type_axis.values = {'Observation'};
        
        SW_xPlt.axis(end + 1) = observed_data_type_axis;
        
        SW_ShuffleMerged = SW_xPlt.merge(SW_Shuffle);
        
        switch norm_struct.how{2}
            
            case ''
        
                SW_xPlt = norm_axis_by_value(SW_ShuffleMerged, 'Data_Type', 'Surrogate');
                
            case 'subtract'
        
                SW_xPlt = norm_axis_by_value(SW_ShuffleMerged, 'Data_Type', 'Surrogate', 'subtract');
                
            case 'zscore'
        
                SW_xPlt = norm_axis_by_value(SW_ShuffleMerged, 'Data_Type', 'Surrogate', 'subtract');
        
                SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
                
                collapse_struct = struct('function', @nanstd, 'varargin', []);
                
                SW_Shuffle = PD_xPlt_prep_for_norm(SW_Shuffle_struct, SW_xPlt, function_name, norm_pd_names, collapse_struct);
        
                SW_Shuffle.axis(end + 1) = surrogate_data_type_axis;
                
                SW_ShuffleMerged = SW_xPlt.merge(SW_Shuffle);
                
                SW_xPlt = norm_axis_by_value(SW_ShuffleMerged, 'Data_Type', 'Surrogate', 'divide');
                
        end
        
        SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
        
        SW_xPlt.getaxisinfo
            
end

end