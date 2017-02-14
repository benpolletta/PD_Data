function PD_sliding_window_pre_post(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin)
    
close('all')

pd_names = {'pre', 'post'};

for s = 1:length(subjects_mat_cell)
        
    subjects_mat_name = subjects_mat_cell{s};
    
    load(subjects_mat_name)
    
    no_folders = length(folders);
    
    subjects_mat_prefix = subjects_mat_name(1:(end - length('_subjects.mat')));
    
    parfor fo = 1:no_folders
        
        subjects_struct = load(subjects_mat_name);
        
        folder = subjects_struct.folders{fo};
        
        prefix = subjects_struct.prefixes{fo};
        
        subj_name = [folder, '/', prefix];
        
        data = load([subj_name, get_file_suffix(data_labels_struct, data_labels_struct.data_field)]);
        
        max_beta_density = load([subjects_mat_prefix, get_file_suffix(data_labels_struct, 'high_density_periods')]);
        
        striatal_channel = find(strcmp(subjects_struct.chan_labels, 'Striatum'));
        
        if isfield(data, data_labels_struct.data_field)
            
            data = data.(data_labels_struct.data_field);
            
        else
            
            display(sprintf('Data not found for %s.', folder))
            
        end
        
        t = (1:size(data, 1))/data_labels_struct.sampling_freq{1} - subjects_struct.basetimes(fo); % t = t/60;
        
        for pd = 1:2
            
            bp_max_start = max_beta_density.All_bp_max_start(fo, striatal_channel, data_labels_struct.band_index, pd);
            
            bp_max_end = max_beta_density.All_bp_max_end(fo, striatal_channel, data_labels_struct.band_index, pd);
            
            max_beta_density_index = t >= t(bp_max_start) & t <= t(bp_max_end);
            
            data_selected = data(max_beta_density_index, :);
            
            [~, ~] = sliding_window_analysis_multichannel(function_handle, data_selected,...
                data_labels_struct.sampling_freq, sliding_window_cell, 1, 1, [subj_name, '_', pd_names{pd},...
                '_band', num2str(data_labels_struct.band_index)], varargin{:});
            
        end
        
    end
    
end

end