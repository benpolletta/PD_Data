function beta_blocks_consolidate(subject_mat, time_window, percent, freqs, no_cycles, bands)

PD_struct = PD_initialize(subject_mat);

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

for fo = 1:PD_struct.no_folders
    
    folder = PD_struct.folders{fo};
    
    prefix = PD_struct.prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    BP_high_consolidated = struct();
    
    for n = 1:PD_struct.no_norms
        
        for ty = 1:PD_struct.no_types
    
            BP_high_data = load([subj_name, BP_suffix, '_2sd_BP', PD_struct.norms{n}, '_high.mat'], ['BP_high', PD_struct.high_type{ty}]);
            
            BP_high_data = getfield(BP_high_data, ['BP_high', PD_struct.high_type{ty}]);
            
            BP_high_data(:, :, 3) = BP_high_data(:, :, 1) | BP_high_data(:, :, 2);
            
            BP_high_data(:, :, 4) = BP_high_data(:, :, 1) & BP_high_data(:, :, 2);
            
            BP_high_cons = nan(size(BP_high_data));
            
            for ch = 1:size(BP_high_data, 3)
                
                for b = 1:PD_struct.no_bands
                    
                    BP_high_conv = conv(BP_high_data(:, b, ch), ones(time_window, 1)/time_window, 'same');
                    
                    BP_high_cons(:, b, ch) = BP_high_conv >= percent;
                    
                end
            
            end
            
            BP_high_consolidated = setfield(BP_high_consolidated, ['BP', PD_struct.norms{n}, '_high', PD_struct.high_type{ty}], BP_high_cons);
            
        end
        
    end
    
    save([subj_name, BP_suffix, '_2sd_BP_high_', num2str(time_window/PD_struct.sampling_freq), 's_', num2str(percent), 'pct_consolidated.mat'], '-v7.3', 'BP_high_consolidated')
    
end

