function PMTM = PD_PMTM_pre_post_load(subject_mat_cell, band_index, filename)
    
close('all')

period_names = {'pre', 'post'}; no_periods = length(period_names);

sampling_freq = 500;

for s = 1:length(subject_mat_cell)
    
    load(subject_mat_cell{s})
    
    subject_mat_folders(s) = length(folders);
    
end

no_folders = sum(subject_mat_folders);

channel_labels = {'Striatum'; 'Motor Ctx.'}; no_channels = length(channel_labels);

no_directions = length(channel_labels);

PMTM = nan(150*sampling_freq + 1, no_folders, no_periods, no_channels);

folder_index = 0;

for s = 1:length(subject_mat_cell)
    
    subjects_struct = load(subject_mat_cell{s});
    
    chan_order = [1 2];
    
    if ~strcmp(chan_labels{1}, 'Striatum');
        
        chan_order = fliplr(chan_order);
        
    end
    
    for fo = 1:length(subjects_struct.folders)
        
        folder_index = folder_index + 1;
        
        folder = subjects_struct.folders{fo};
        
        prefix = subjects_struct.prefixes{fo};
        
        subj_name = [folder, '/', prefix];
        
        for period = 1:no_periods
            
            analysis_name = make_sliding_window_analysis_name([subj_name, '_', period_names{period},...
                '_band', num2str(band_index)],'pmtm',{[150 150],[1 1]},2);
            
            load(analysis_name)
            
            for channel = 1:no_channels
            
                PMTM(:, folder_index, period, channel) = sw(:, chan_order(channel));
                
            end
            
        end
        
    end
    
end

save(make_sliding_window_analysis_name([filename, '_pre_post_band', num2str(band_index)],...
    'pmtm',{[150 150],[1 1]},2), 'PMTM')