function run_phase_analysis(subject_mat, peak_suffix, time_window, percent, f_tol, freqs, no_cycles, bands, band_index, plot_opt, epoch_secs_for_plot)

% time_window is the length of the window in which you look for at minimum 
%   percent of data to be over 2 s.d. for beta power (in samples).
% percent is proportion of data that should be over 2 s.d. for beta power
%   in given time window (of length time_window) (given as a proportion).

beta_blocks_consolidate(subject_mat, peak_suffix, time_window, percent, freqs, no_cycles, bands)

if plot_opt == 1

    beta_blocks_consolidated_plot(subject_mat, peak_suffix, epoch_secs_for_plot, time_window, percent, band_index, freqs, no_cycles, bands)

end
    
beta_blocks_consolidated_freqs(subject_mat, peak_suffix, time_window, percent, '', band_index, freqs, no_cycles, bands)

beta_blocks_consolidated_phase_analysis(subject_mat, peak_suffix, time_window, percent, '', band_index, f_tol, freqs, no_cycles, bands)