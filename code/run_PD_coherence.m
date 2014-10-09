function run_PD_coherence(subject_mat, outlier_lim, sd_lim, win_size, smooth_size, tbw)

% PD_beta_epochs_rel_infusion(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)
% 
% close('all')

% PD_beta_epochs_coh_mtm(subject_mat, outlier_lim, sd_lim, win_size, smooth_size, tbw)
% 
% close('all')

PD_beta_epochs_coh_mtm_collect(subject_mat, outlier_lim, sd_lim, win_size, smooth_size, tbw)

close('all')

PD_beta_epochs_coh_mtm_plot_group(subject_mat, outlier_lim, sd_lim, win_size, smooth_size, tbw)

close('all')

% PD_beta_epochs_coh_mtm_plot_individual(subject_mat, outlier_lim, sd_lim, win_size, smooth_size, tbw)
% 
% close('all')
