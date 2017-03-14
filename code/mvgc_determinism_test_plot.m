function mvgc_determinism_test_plot(subjects_mat_prefix)

subjects_mat_name = [subjects_mat_prefix, '_subjects.mat'];

% subjects_mat_prefix = subjects_mat_name(1:(end - length('_subjects.mat')));

subjects_struct = load(subjects_mat_name);

striatal_channel = find(strcmp(subjects_struct.chan_labels, 'Striatum'));

figure();

for fo = 1:length(subjects_struct.folders)
    
    folder = subjects_struct.folders{fo};
    
    prefix = subjects_struct.prefixes{fo};
    
    subj_name = [folder, '/', prefix];
    
    load([subj_name, '_mvgc_determinism_test.mat'])
    
    f = permute(f, [1 4 2 3 5]);
    
    frequencies = 250*(1:size(f, 1))/size(f,1);
    
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

save_as_pdf(gcf, [subjects_mat_prefix, '_mvgc_determinism_test'])