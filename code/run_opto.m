function run_opto(directory, subject_mat, freqs, no_cycles, bands)

present_dir = pwd;

if ~strcmp(present_dir((end - length(directory) + 1):end), directory)
    
    cd (directory)
    
end

% PD_laser_artifacts_wav(subject_mat)

PD_beta_blocks_rel_infusion(subject_mat, 2, [], freqs, no_cycles, bands)

pd_handles = {'', '_power', '_power'};

norm_handles = {'', '', '_pct'};

for i = 2:length(pd_handles)

    PD_beta_blocks_rel_infusion_laser_stats(subject_mat, pd_handles{i}, norm_handles{i}, 10, freqs, no_cycles, bands)
    
    PD_beta_blocks_rel_infusion_laser_plots(subject_mat, pd_handles{i}, 10, freqs, no_cycles, bands)
    
    PD_beta_blocks_rel_infusion_laser_plots_individual(subject_mat, pd_handles{i}, norm_handles{i}, 10, freqs, no_cycles, bands)
    
end