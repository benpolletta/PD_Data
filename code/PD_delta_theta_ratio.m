function delta_theta = PD_delta_theta_ratio

sb = load('seven_bands');

data_labels_struct = init_data_labels(sb.freqs, sb.no_cycles, sb.bands);

matnames = {'st_m1', 'st_m1_ali', 'st_m1_ali2'};

norm = '';

folder_index = 1;

for m = 1:length(matnames)
    
    sm = load([matnames{m}, '_subjects.mat']);
    
    for f = 1:length(sm.folders)
        
        clear BP
        
        BP = get_BP([sm.folders{f}, '/', sm.prefixes{f}], data_labels_struct.peak_suffix, [], sm.outlier_lims(f), norm, sb.freqs, sb.no_cycles, sb.bands);
        
        time = (1:size(BP, 1))/500 - sm.basetimes(f);
        
        delta_theta(folder_index, :) = nanmean(BP(time <= 0, 1, :)./BP(time <= 0, 2, :));
        
        folder_index = folder_index + 1;
    
    end
    
end

save('PD_delta_theta_ratio.mat', 'delta_theta')

end