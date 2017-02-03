function increase_summary = run_power_post_vs_pre(subject_mat_cell, save_name, channel_labels, peak_suffix, freqs, no_cycles, bands, band_indices, window_length)

[freqs, no_cycles, bands, ~] = init_freqs(freqs, no_cycles, bands);

% for s = 1:length(subject_mat_cell)
%     
%     power_post_vs_pre_stats(subject_mat_cell{s}, peak_suffix, 150, '_by_STR', '_pct', freqs, no_cycles, bands, band_indices, window_length)
% 
% end
    
increase_summary = power_post_vs_pre_summary(subject_mat_cell, save_name,...
    channel_labels, peak_suffix, '_pct', freqs, no_cycles, bands, band_indices, window_length);

power_post_vs_pre_timeseries_data(subject_mat_cell, save_name, peak_suffix,...
    '_pct', freqs, no_cycles, bands, band_indices, window_length)

power_post_vs_pre_timeseries_plot(save_name, channel_labels, peak_suffix, freqs, no_cycles, bands, band_indices, window_length)
