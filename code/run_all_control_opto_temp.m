function run_all_control_opto(peak_suffix, freqs, no_cycles, bands, band_index)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    no_bands = 6;
    
else
    
    no_bands = size(bands, 1);
    
end

subject_mat_prefixes = {'st_m1_ali_control_opto', 'st_m1_ali2_control_opto', 'st_m1_ali3_control_opto'};

no_mats = length(subject_mat_prefixes);

for m = 2:no_mats

    run_opto('Control_Opto', [subject_mat_prefixes{m}, '_subjects.mat'], peak_suffix, freqs, no_cycles, bands)

end

measures = {'', '_power'}; norms = {'', '_pct'};

for m = 1:length(measures)
    
    collect_control_opto_power_density(peak_suffix, measures{m}, norms{m}, freqs, no_cycles, bands)
    
    PD_beta_blocks_rel_infusion_pre_post_plot('CONTROL_OPTO_subjects.mat', peak_suffix, [], ['_10trials', norms{m}, measures{m}], 't-test', freqs, no_cycles, bands)

end
   
collect_control_opto_spectrum(peak_suffix, '_pct', no_bands, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion_pre_post_spectrum_plot('CONTROL_OPTO_subjects.mat', peak_suffix, [], '_10trials', '_pct', no_bands, no_bands, freqs, no_cycles, bands)

make_opto_spectra_figures('CONTROL_OPTO', peak_suffix)

for m = 1:no_mats

    for ftol = 2:3
        
        run_phase_analysis([subject_mat_prefixes{m}, '_subjects.mat'], peak_suffix, 200, .5, ftol, freqs, no_cycles, bands, band_index, 1, [])
        
    end

end

collect_control_opto_freqs(peak_suffix, 200, .5, '', band_index, freqs, no_cycles, bands)

for ftol = 2:3
    
    beta_blocks_consolidated_phase_analysis('CONTROL_OPTO_subjects.mat', peak_suffix, 200, .5, '', band_index, ftol, freqs, no_cycles, bands)
    
end