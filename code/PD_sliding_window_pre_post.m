function PD_sliding_window_pre_post(function_handle, sampling_freq, sliding_window_cell, subjects_mat, data_labels_struct, band_index, varargin)
    
close('all')

load(subjects_mat)

subjects_mat_name = subjects_mat(1:(end - length('_subjects.mat')));

no_folders = length(folders);

pd_names = {'pre', 'post'};

parfor fo = 1:no_folders
    
    subjects_struct = load(subjects_mat);
    
    folder = subjects_struct.folders{fo};
    
    prefix = subjects_struct.prefixes{fo};
    
    subj_name = [folder, '/', prefix];
    
    data = load([subj_name, get_file_suffix(data_labels_struct, data_labels_struct.data_field)]);

    max_beta_density = load([subjects_mat_name, get_file_suffix(data_labels_struct, 'high_density_periods')]);
    
    subjects_struct = load(subjects_mat);
    
    folder = subjects_struct.folders{fo};
    
    prefix = subjects_struct.prefixes{fo};
    
    striatal_channel = find(strcmp(subjects_struct.chan_labels, 'Striatum'));
    
    subj_name = [folder,'/',prefix];
    
    if isfield(data, data_labels_struct.data_field)
        
        data = data.(data_labels_struct.data_field);
        
    else
        
        display(sprintf('Data not found for %s.', folder))
        
    end
    
    t = (1:size(data, 1))/sampling_freq{1} - subjects_struct.basetimes(fo); % t = t/60;
    
    for pd = 1:2
        
        bp_max_start = max_beta_density.All_bp_max_start(fo, striatal_channel, band_index, 2);
        
        bp_max_end = max_beta_density.All_bp_max_end(fo, striatal_channel, band_index, 2);
        
        max_beta_density_index = t >= t(bp_max_start) & t <= t(bp_max_end);
        
        data_selected = data(max_beta_density_index, :);
        
        [~, ~] = sliding_window_analysis_multichannel(function_handle, data_selected,...
            sampling_freq, sliding_window_cell, 1, 1, [subj_name, pd_names{pd}, '_band', num2str(band_index)], varargin{:});
        
    end
    
end

end