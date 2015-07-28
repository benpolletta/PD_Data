function plot_carb_spectrum(norm, band_index_for_time, band_index_for_display, freqs, no_cycles, bands)

matnames = {'st_m1', 'st_stn', 'st_m1_6OHDA'};

period_labels = {'_by_STR', '_by_STR', ''};

matlabpool open 3

parfor i = 1:3

    PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual([matnames{i}, '_subjects.mat'],...
        150, period_labels{i}, norm, band_index_for_time, band_index_for_display, freqs, no_cycles, bands)
    
end

matlabpool close

PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual('STR_subjects.mat',...
    150, '', norm, band_index_for_time, band_index_for_display, freqs, no_cycles, bands)

PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual('CARB_6OHDA_subjects.mat',...
    150, '', norm, band_index_for_time, band_index_for_display, freqs, no_cycles, bands)