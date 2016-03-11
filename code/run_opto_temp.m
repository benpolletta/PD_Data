function run_opto(directory, subject_mat, peak_suffix, freqs, no_cycles, bands)

present_dir = pwd;

if ~strcmp(present_dir((end - length(directory) + 1):end), directory)
    
    cd (directory)
    
end

% PD_pct_fix(subject_mat, peak_suffix, freqs, no_cycles, bands)

% PD_beta_blocks_rel_infusion(subject_mat, peak_suffix, 2, freqs, no_cycles, bands)

pd_handles = {'', '_power', '_power'};

norm_handles = {'', '', '_pct'};

% for i = 1:length(pd_handles)
% 
%      PD_beta_blocks_rel_infusion_laser_stats(subject_mat, peak_suffix, pd_handles{i}, norm_handles{i}, 10, freqs, no_cycles, bands)
%      
%      PD_beta_blocks_rel_infusion_laser_plots(subject_mat, peak_suffix, pd_handles{i}, 10, freqs, no_cycles, bands)
%      
%      PD_beta_blocks_rel_infusion_laser_plots_individual(subject_mat, peak_suffix, pd_handles{i}, norm_handles{i}, 10, freqs, no_cycles, bands)
%     
% end

for i = 1:2
   
    PD_beta_blocks_rel_infusion_laser_spectrum(subject_mat, peak_suffix, norm_handles{i + 1}, 10, 6, freqs, no_cycles, bands)
    
    PD_beta_blocks_rel_infusion_pre_post_spectrum_plot(subject_mat, peak_suffix, [], '_10trials', norm_handles{i + 1}, 6, 6, freqs, no_cycles, bands);
    
    PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual(subject_mat, peak_suffix, [], '_10trials', norm_handles{i + 1}, 6, 6, freqs, no_cycles, bands);
    
end

PD_rel_infusion_plot_spectrogram(subject_mat, peak_suffix, 2, 20)