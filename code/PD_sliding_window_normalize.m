function SW_xPlt = PD_sliding_window_normalize(SW_xPlt, function_name, sliding_window_cell, data_labels_struct, filename, pd_names, norm_struct, norm_pd_names, varargin)

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
        if ~strcmp(function_name, 'PAC')
            
            frequencies = SW_xPlt.meta.matrix_dim_1.values;
            
        else
            
            error('Frequency normalization has not been coded up for PAC data.')
            
        end
        
        freq_mult = diag(frequencies.^(2/3));
        
        SW_xPlt.data = cellfun(@(x) freq_mult*x, SW_xPlt.data, 'UniformOutput', 0);
        
    case 'baseline'
    %% Normalization by mean baseline value.
    
        % Getting name of file from which to load baseline data.
        if isempty(norm_pd_names)
            
            baseline_pd_label = '_baseline';
            
            norm_pd_names = {'baseline'};
            
        else
            
            baseline_pd_label = '';
            
            for pd = 1:length(norm_pd_names)
                
                baseline_pd_label = [baseline_pd_label, '_', norm_pd_names{pd}];
                
            end
            
        end

        SW_Baseline_struct = load([make_sliding_window_analysis_name([filename, baseline_pd_label, '_band',... % '_baseline_band',... % 
            num2str(data_labels_struct.band_index), make_label('bands', legnth(data_labels_struct.bands), 7, 'back')], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        collapse_struct = struct('function', @nanmean, 'varargin', {});
        
        SW_Baseline = PD_xPlt_prep_for_norm(SW_Baseline_struct, SW_xPlt, function_name, norm_pd_names, collapse_struct);
        
        SW_BaselineMerged = SW_xPlt.merge(SW_Baseline);
        
        SW_BaselineMerged.getaxisinfo
        
        switch norm_struct.how
            
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
        
        % SW_xPlt.data = cellfun(@(x, y) 100*(x./y - 1), SW_xPlt.data, SW_Baseline.data, 'UniformOutput', 0);
        
    case 'shuffle'
    %% Normalization by shuffled surrogate data.
        
        if isempty(norm_pd_names)
            
            shuffle_pd_label = '';
            
            for pd = 1:length(pd_names)
                
                norm_pd_names{pd} = [pd_names{pd}, '_shuffles'];
                
                shuffle_pd_label = [shuffle_pd_label, '_', norm_pd_names{pd}];
                
            end
            
        else
            
            shuffle_pd_label = '';
            
            for pd = 1:length(norm_pd_names)
                
                shuffle_pd_label = [shuffle_pd_label, '_', norm_pd_names{pd}];
                
            end
            
        end

        SW_Shuffle_struct = load([make_sliding_window_analysis_name([filename, shuffle_pd_label, '_band',...
            num2str(data_labels_struct.band_index), make_label('bands', legnth(data_labels_struct.bands), 7, 'back')], function_name,...
            window_time_cell, 2, varargin{:}), '_xPlt.mat']);
        
        collapse_struct = struct('function', @nanmean, 'varargin', {});
        
        SW_Shuffle = PD_xPlt_prep_for_norm(SW_Shuffle_struct, SW_xPlt, function_name, norm_pd_names, collapse_struct);
        
        surrogate_data_type_axis = nDDictAxis;
        
        surrogate_data_type_axis.name = 'Data_Type'; surrogate_data_type_axis.values = {'Surrogate'};
        
        SW_Shuffle.axis(end + 1) = surrogate_data_type_axis;
        
        observed_data_type_axis = nDDictAxis;
        
        observed_data_type_axis.name = 'Data_Type'; observed_data_type_axis.values = {'Observation'};
        
        SW_xPlt.axis(end + 1) = observed_data_type_axis;
        
        SW_ShuffleMerged = SW_xPlt.merge(SW_Shuffle);
        
        switch norm_struct.how
            
            case ''
        
                SW_xPlt = norm_axis_by_value(SW_ShuffleMerged, 'Data_Type', 'Surrogate');
                
            case 'subtract'
        
                SW_xPlt = norm_axis_by_value(SW_ShuffleMerged, 'Data_Type', 'Surrogate', 'subtract');
                
            case 'zscore'
        
                SW_xPlt = norm_axis_by_value(SW_ShuffleMerged, 'Data_Type', 'Surrogate', 'subtract');
        
                SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
                
                collapse_struct = struct('function', @nanstd, 'varargin', []);
                
                SW_Shuffle_struct.SW_xPlt.axis(end + 1) = surrogate_data_type_axis;
                
                SW_Shuffle = PD_xPlt_prep_for_norm(SW_Shuffle_struct, SW_xPlt, function_name, norm_pd_names, collapse_struct);
                
                SW_ShuffleMerged = SW_xPlt.merge(SW_Shuffle);
                
                SW_xPlt = norm_axis_by_value(SW_ShuffleMerged, 'Data_Type', 'Surrogate', 'divide');
                
        end
        
        SW_xPlt = SW_xPlt.axissubset('Data_Type', 'Observation');
        
        SW_xPlt.getaxisinfo
            
end

end