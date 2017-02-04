function PD_sliding_GC_baseline_post(window_length, subjects_mat, peak_suffix, data_suffix, epoch_secs, pd_handle, freqs, no_cycles, bands, band_index)

[~, ~, ~, BP_suffix] = init_freqs(freqs, no_cycles, bands);
    
close('all')

load(subjects_mat)

subjects_mat_name = subjects_mat(1:(end - length('_subjects.mat')));

no_folders = length(folders);

calc_names = {'baseline', 'post'};

sampling_freq = 500;

parfor fo = 1:no_folders

    max_beta_density = load([subjects_mat_name, BP_suffix, peak_suffix,...
        '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, '.mat']);
    
    subjects_struct = load(subjects_mat);
    
    folder = subjects_struct.folders{fo};
    
    prefix = subjects_struct.prefixes{fo};
    
    striatal_channel = find(strcmp(subjects_struct.chan_labels, 'Striatum'));
    
    subj_name = [folder,'/',prefix];
    
    data = load([subj_name, '_all_channel_data_dec', data_suffix, '.mat']);
    
    if isfield(data, 'PD_dec')
        
        data = data.PD_dec;
        
    elseif isfield(data, 'data_subtracted')
        
        data = data.data_subtracted;
        
    else
        
        display(sprintf('Data not found for %s.', folder))
        
    end
    
    t = (1:size(data, 1))/sampling_freq - subjects_struct.basetimes(fo); % t = t/60;
        
    baseline_index = t < 0; 
    
    baseline_index = logical(baseline_index);
    
    bp_max_start = max_beta_density.All_bp_max_start(fo, striatal_channel, band_index, 2);
    
    bp_max_end = max_beta_density.All_bp_max_end(fo, striatal_channel, band_index, 2);
    
    max_beta_density_index = t >= t(bp_max_start) & t <= t(bp_max_end);
    
    indices_calculated = [baseline_index' max_beta_density_index'];
    
    no_calcs = size(indices_calculated, 2);
    
    for calc = 1:no_calcs
        
        data_selected = data(indices_calculated(:, calc), :);
        
        sliding_window_analysis_multichannel(@mvgc_analysis, data_selected,...
            sampling_freq, window_length, window_length, 1, 1, [subj_name, calc_names{calc}], [], [], []);
        
    end
    
end

end