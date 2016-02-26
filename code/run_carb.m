function run_carb(subject_mat, peak_suffix, freqs, no_cycles, bands)

PD_pct_fix(subject_mat, peak_suffix, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion(subject_mat, peak_suffix, 2, freqs, no_cycles, bands);

PD_beta_blocks_rel_infusion_pre_post_pds(subject_mat, peak_suffix, 150, freqs, no_cycles, bands);

power_handles = {'', '_pct'};

for i = 1:length(power_handles)
    
    PD_beta_blocks_rel_infusion_pre_post_power(subject_mat, peak_suffix, 150, '_by_STR', power_handles{i}, freqs, no_cycles, bands);
   
    PD_beta_blocks_rel_infusion_pre_post_spectrum(subject_mat, peak_suffix, 150, '_by_STR', power_handles{i}, 3, freqs, no_cycles, bands);

end

pd_handles = {'_by_STR', '_by_STR_power', '_by_STR_pct_power'};

for i = 1:length(pd_handles)
    
    PD_beta_blocks_rel_infusion_pre_post_plot(subject_mat, peak_suffix, 150, ['_by_STR', pd_handles{i}], 'ranksum', freqs, no_cycles, bands);
    
    PD_beta_blocks_rel_infusion_pre_post_plot_individual(subject_mat, peak_suffix, 150, 'ranksum', ['_by_STR', pd_handles{i}], freqs, no_cycles, bands);
   
end

pd_handles = {'_by_STR', '_by_STR_pct'};

for i = 1:length(pd_handles)
    
    PD_beta_blocks_rel_infusion_pre_post_spectrum_plot(subject_mat, peak_suffix, 150, pd_handles{i}, power_handles{j}, 3, 6, freqs, no_cycles, bands);
    
    PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual(subject_mat, peak_suffix, 150, pd_handles{i}, power_handles{j}, 3, 6, freqs, no_cycles, bands);
    
end

PD_rel_infusion_plot_spectrogram(subject_mat, peak_suffix, 2, 20)