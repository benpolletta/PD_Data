function run_opto(directory, subject_mat, peak_suffix, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    no_bands = 6;
    
else
    
    no_bands = size(bands, 1);
    
end

present_dir = pwd;

if ~strcmp(present_dir((end - length(directory) + 1):end), directory)
    
    cd (directory)
    
end

PD_pct_fix(subject_mat, peak_suffix, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion(subject_mat, peak_suffix, 2, freqs, no_cycles, bands)

pd_handles = {'', '_power', '_power'};

norm_handles = {'', '', '_pct'};

for i = 1:length(pd_handles)

    PD_beta_blocks_rel_infusion_laser_stats(subject_mat, peak_suffix, pd_handles{i}, norm_handles{i}, 10, freqs, no_cycles, bands)
    
    PD_beta_blocks_rel_infusion_laser_plots(subject_mat, peak_suffix, pd_handles{i}, 10, freqs, no_cycles, bands)
    
    PD_beta_blocks_rel_infusion_laser_plots_individual(subject_mat, peak_suffix, pd_handles{i}, norm_handles{i}, 10, freqs, no_cycles, bands)
    
end

for i = 1:2
   
    PD_beta_blocks_rel_infusion_laser_spectrum(subject_mat, peak_suffix, norm_handles{i + 1}, 10, no_bands, freqs, no_cycles, bands)
    
    PD_beta_blocks_rel_infusion_pre_post_spectrum_plot(subject_mat, peak_suffix, [], '_10trials', norm_handles{i + 1}, no_bands, no_bands, freqs, no_cycles, bands);
    
    PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual(subject_mat, peak_suffix, [], '_10trials', norm_handles{i + 1}, no_bands, no_bands, freqs, no_cycles, bands);
    
end

PD_beta_blocks_rel_infusion_laser_PLV(subject_mat, peak_suffix, 10, band_index, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion_pre_post_PLV_plot(subject_mat, peak_suffix, [], '_10trials', no_bands, no_bands, freqs, no_cycles, bands)

PD_rel_infusion_plot_spectrogram(subject_mat, peak_suffix, 2, 20)