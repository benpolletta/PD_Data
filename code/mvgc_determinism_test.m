function mvgc_determinism_test(subjects_mat_prefix, number_iterations)

freqs = load('seven_bands.mat', 'freqs'); freqs = freqs.freqs;

no_cycles = load('seven_bands.mat', 'no_cycles'); no_cycles = no_cycles.no_cycles;

bands = load('seven_bands.mat', 'bands'); bands = bands.bands;

data_labels_struct = init_data_labels(freqs, no_cycles, bands, 'data_field', 'data_subtracted');

subjects_mat_name = [subjects_mat_prefix, '_subjects.mat'];

% subjects_mat_prefix = subjects_mat_name(1:(end - length('_subjects.mat')));

subjects_struct = load(subjects_mat_name);

striatal_channel = find(strcmp(subjects_struct.chan_labels, 'Striatum'));

for fo = 1:length(subjects_struct.folders)
    
    folder = subjects_struct.folders{fo};
    
    prefix = subjects_struct.prefixes{fo};
    
    subj_name = [folder, '/', prefix];
    
    data = load([subj_name, get_file_suffix(data_labels_struct, data_labels_struct.data_field)]);
    
    max_beta_density = load([subjects_mat_prefix, get_file_suffix(data_labels_struct, 'high_density_periods')]);
    
    if isfield(data, data_labels_struct.data_field)
        
        data = data.(data_labels_struct.data_field);
        
        display(sprintf('Evaluating %s for %s.', func2str(@mvgc_analysis), folder))
        
    else
        
        display(sprintf('Data not found for %s.', folder))
        
    end
    
    t = (1:size(data, 1))/data_labels_struct.sampling_freq{1} - subjects_struct.basetimes(fo); % t = t/60;
    
    f = nan(75002, 2, 2, number_iterations, 2);
    
    moAIC = nan(number_iterations, 2);
    
    info = cell(number_iterations, 2);
    
    for pd = 1:2
        
        bp_max_start = max_beta_density.All_bp_max_start(fo, striatal_channel, data_labels_struct.band_index, pd);
        
        bp_max_end = max_beta_density.All_bp_max_end(fo, striatal_channel, data_labels_struct.band_index, pd);
        
        max_beta_density_index = t >= t(bp_max_start) & t <= t(bp_max_end);
        
        data_selected = data(max_beta_density_index, :);
        
        parfor i = 1:number_iterations
            
            display(sprintf('Iteration %d.', i))
            
            [fi, moAICi, infoi] = mvgc_analysis(data_selected, [], '', 1);
            
            % size(fi), size(moAICi), size(infoi)
            
            f(:, :, :, i, pd) = fi;
            
            moAIC(i, pd) = moAICi;
            
            info{i, pd} = infoi;
            
        end
        
    end
    
    save([subj_name, '_mvgc_determinism_test'], 'f', 'moAIC', 'info')
    
end