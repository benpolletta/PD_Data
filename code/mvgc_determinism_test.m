function mvgc_determinism_test(number_iterations)

freqs = load('seven_bands.mat', 'freqs'); freqs = freqs.freqs;

no_cycles = load('seven_bands.mat', 'no_cycles'); no_cycles = no_cycles.no_cycles;

bands = load('seven_bands.mat', 'bands'); bands = bands.bands;

data_labels_struct = init_data_labels(freqs, no_cycles, bands, 'data_field', 'data_subtracted');

subjects_mat_name = 'st_m1_subjects.mat';

subjects_mat_prefix = 'st_m1';

subjects_struct = load(subjects_mat_name);

striatal_channel = find(strcmp(subjects_struct.chan_labels, 'Striatum'));

f = nan(75002, 2, 2, number_iterations, 2);

moAIC = nan(number_iterations, 2);

info = cell(number_iterations, 2);

figure

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
    
    for pd = 1:2
        
        bp_max_start = max_beta_density.All_bp_max_start(fo, striatal_channel, data_labels_struct.band_index, pd);
        
        bp_max_end = max_beta_density.All_bp_max_end(fo, striatal_channel, data_labels_struct.band_index, pd);
        
        max_beta_density_index = t >= t(bp_max_start) & t <= t(bp_max_end);
        
        data_selected = data(max_beta_density_index, :);
        
        parfor i = 1:number_iterations
            
            display(sprintf('Iteration %d.', i))
            
            [fi, moAICi, infoi] = mvgc_analysis(data_selected, [], '', 1);
            
            size(fi), size(moAICi), size(infoi)
            
            f(:, :, :, i, pd) = fi;
            
            moAIC(i, pd) = moAICi;
            
            info{i, pd} = infoi;
            
        end
        
    end
    
    save([folder, 'mvgc_determinism_test'], 'f', 'moAIC', 'info')
    
    f = permute(f, [1 4 2 3 5]);
    
    frequencies = 500*(1:size(f, 1))/size(f,1);
    
    dir_channels = {1, 2};
    
    for d = 1:2
        
        dir_channels = fliplr(dir_channels);
        
        subplot(length(subjects_struct.folders), 2, (fo - 1)*2 + d)
        
        % subplot_index = (p - 1)*2 + d;
        %
        % subplot(2, 2, subplot_index)
        
        f_mean = squeeze(nanmean(f(:, :, dir_channels{:}, :), 2));
        
        f_std = squeeze(nanstd(f(:, :, dir_channels{:}, :), [], 2));
        
        boundedline(frequencies, f_mean, prep_for_boundedline(f_std))
        
        box off
        
        ylabel(folder)
        
    end
    
end

save_as_pdf(gcf, 'mvgc_determinism_test')