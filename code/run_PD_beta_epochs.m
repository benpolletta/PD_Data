function run_PD_beta_epochs(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

PD_beta_epochs_rel_infusion(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

close('all')

PD_beta_epochs_rel_infusion_roseplot_by_datapoint(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

close('all')

PD_beta_epochs_rel_infusion_roseplot_by_datapoint_group(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

close('all')

PD_beta_epochs_rel_infusion_roseplot_by_datapoint_individual(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

close('all')

PD_beta_epochs_rel_infusion_roseplot_by_dp_ind_avg(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)