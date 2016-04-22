function run_all_carb_opto(peak_suffix, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    no_bands = 6;
    
else
    
    no_bands = size(bands, 1);
    
end

subject_mat_prefixes = {'st_m1_emxarch', 'st_m1_ali_post_carb_opto', 'st_m1_ali2_post_carb_opto'};

no_mats = length(subject_mat_prefixes);

for m = 1:no_mats

    run_opto('Carb_Opto', [subject_mat_prefixes{m}, '_subjects.mat'], peak_suffix, freqs, no_cycles, bands)

end

measures = {'', '_power'}; norms = {'', '_pct'};

for m = 1:length(measures)
    
    collect_opto_power_density(peak_suffix, measures{m}, norms{m}, freqs, no_cycles, bands)
    
    PD_beta_blocks_rel_infusion_pre_post_plot('OPTO_subjects.mat', peak_suffix, [], ['_10trials', norms{m}, measures{m}], 't-test', freqs, no_cycles, bands)

end
   
collect_opto_spectrum(peak_suffix, '_pct', no_bands, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion_pre_post_spectrum_plot('OPTO_subjects.mat', peak_suffix, [], '_10trials', '_pct', no_bands, no_bands, freqs, no_cycles, bands)

make_opto_spectra_figures('OPTO', peak_suffix)
   
collect_opto_PLV(peak_suffix, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion_pre_post_PLV_plot('OPTO_subjects.mat', peak_suffix, [], '_10trials', no_bands, no_bands, freqs, no_cycles, bands)

make_opto_PLV_figures('OPTO', peak_suffix)

% for m = 1:no_mats
%         
%     run_phase_analysis([subject_mat_prefixes{m}, '_subjects.mat'], peak_suffix, 200, .5, ftol, freqs, no_cycles, bands, band_index, 1, [])
% 
% end
% 
% collect_opto_freqs(peak_suffix, 200, .5, '', band_index, freqs, no_cycles, bands)
% 
% for ftol = 2:3
%     
%     beta_blocks_consolidated_phase_analysis('OPTO_subjects.mat', peak_suffix, 200, .5, '', band_index, ftol, freqs, no_cycles, bands)
%     
% end