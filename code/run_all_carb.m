function run_all_carb(peak_suffix, freqs, no_cycles, bands, band_index)

subject_mat_prefixes = {'st_m1', 'st_m1_ali', 'st_m1_ali2'};

no_mats = length(subject_mat_prefixes);

for m = 1:no_mats

    run_carb([subject_mat_prefixes{m}, '_subjects.mat'], peak_suffix, freqs, no_cycles, bands)

end

measures = {'', '_power'}; norms = {'', '_pct'};

for m = 1:length(measures)
    
    collect_striatal_motor_power_density(peak_suffix, measures{m}, norms{m}, freqs, no_cycles, bands)
    
    PD_beta_blocks_rel_infusion_pre_post_plot('STR_M1_subjects.mat', peak_suffix, 150, [norms{m}, measures{m}], 't-test', freqs, no_cycles, bands)

end
   
collect_striatal_w_motor_spectrum(peak_suffix, 150, '_pct', 3, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion_pre_post_spectrum_plot('STR_w_M1_subjects.mat', peak_suffix, 150, '_power', '_pct', band_index, 6, freqs, no_cycles, bands)

collect_motor_spectrum(peak_suffix, 150, '_pct', 3, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion_pre_post_spectrum_plot('M1_subjects.mat', peak_suffix, 150, '_power', '_pct', band_index, 6, freqs, no_cycles, bands)

for m = 1:no_mats

    for ftol = 2:3
        
        run_phase_analysis([subject_mat_prefixes{m}, '_subjects.mat'], peak_suffix, 200, .5, ftol, freqs, no_cycles, bands, band_index, 1, 150)
        
    end

end

collect_striatal_motor_freqs(peak_suffix, 200, .5, '', band_index, freqs, no_cycles, bands)