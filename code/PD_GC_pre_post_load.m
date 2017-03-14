function GC = PD_GC_pre_post_load(subject_mat_cell, band_index, filename)
    
close('all')

pd_names = {'pre', 'post'}; no_periods = length(pd_names);

sampling_freq = 500;

for s = 1:length(subject_mat_cell)
    
    load(subject_mat_cell{s})
    
    subject_mat_folders(s) = length(folders);
    
end

no_folders = sum(subject_mat_folders);

direction_labels = {'Str. -> M1'; 'M1 -> Str.'};

no_directions = length(direction_labels);

GC = nan(150*sampling_freq + 1, no_folders, no_periods, no_directions);

folder_index = 0;

for s = 1:length(subject_mat_cell)
    
    subjects_struct = load(subject_mat_cell{s});
    
    chan_order = [1 2];
    
    if ~strcmp(subjects_struct.chan_labels{1}, 'Striatum');
        
        chan_order = fliplr(chan_order);
        
    end
    
    for fo = 1:length(subjects_struct.folders)
        
        folder_index = folder_index + 1;
        
        folder = subjects_struct.folders{fo}
        
        prefix = subjects_struct.prefixes{fo};
        
        subj_name = [folder, '/', prefix];
        
        for pd = 1:2
            
            analysis_name = make_sliding_window_analysis_name([subj_name, '_', pd_names{pd},...
                '_band', num2str(band_index)], 'mvgc_analysis', {[150 150], [2 2]}, 2); % , [], [], []);
            
            load(analysis_name)
            
            pd_names{pd}
            
            for direction = 1:2
                
                dir_indices = fliplr(chan_order);
                
                if direction == 2, dir_indices = fliplr(dir_indices); end
                
                dir_indices
            
                GC(:, folder_index, pd, direction) = sw(:, dir_indices(1), dir_indices(2));
                
            end
            
        end
        
    end
    
end

save(make_sliding_window_analysis_name([filename, '_pre_post_band', num2str(band_index)],...
    'mvgc_analysis', {[150 150], [2 2]}, 2), 'GC')