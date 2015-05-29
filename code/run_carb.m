function run_carb(subject_mat, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion(subject_mat, 2, [], freqs, no_cycles, bands);

PD_beta_blocks_rel_infusion_pre_post_pds(subject_mat, 150, freqs, no_cycles, bands);

power_handles = {'', '_pct'};

for i = 1:length(power_handles)
    
    PD_beta_blocks_rel_infusion_pre_post_power(subject_mat, 150, '_by_STR', power_handles{i}, freqs, no_cycles, bands);
    
end

pd_handles = {'_by_STR', '_by_STR_power', '_by_STR_pct_power'};

for i = 1:length(pd_handles)
    
    PD_beta_blocks_rel_infusion_pre_post_plot(subject_mat, 150, pd_handles{i}, freqs, no_cycles, bands);
    
    PD_beta_blocks_rel_infusion_pre_post_plot_individual(subject_mat, 150, 'ranksum', pd_handles{i}, freqs, no_cycles, bands);
    
end
