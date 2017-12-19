function run_carb_groups(peak_suffix, time_window, freqs, no_cycles, bands, band_indices)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    no_bands = 6;
    
else
    
    no_bands = size(bands, 1);
    
end

load('missing_2')

load(['M1_groups', make_label('win', time_window, [])])

groups = {missing_2, M1_increased, M1_not_increased};
        
peak_suffix = [peak_suffix, make_label('win', time_window, [])];

for band_index = band_indices
    
    for g = 1:length(groups)
        
        PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual('STR_w_M1_subjects.mat', peak_suffix, 150, '', '_pct', band_index, no_bands, freqs, no_cycles, bands, groups{g}, [])
        
        PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual('M1_subjects.mat', peak_suffix, 150, '', '_pct', band_index, no_bands, freqs, no_cycles, bands, groups{g}, [])
        
        PD_beta_blocks_rel_infusion_pre_post_PLV_plot_individual('STR_M1_subjects.mat', peak_suffix, 150, '', band_index, no_bands, freqs, no_cycles, bands, groups{g})
        
    end
    
end