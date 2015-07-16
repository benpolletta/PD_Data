function collect_carb_spectrum(norm, band_index, freqs, no_cycles, bands)

matnames = {'st_m1', 'st_stn', 'st_m1_6OHDA'};

matlabpool open 3

parfor i = 1:3

    PD_beta_blocks_rel_infusion_pre_post_spectrum([matnames{i}, '_subjects.mat'],...
        150, '_by_STR', norm, band_index, freqs, no_cycles, bands)
    
end

collect_striatal_spectrum(150, norm, band_index, freqs, no_cycles, bands)

combine_carb_6OHDA_spectrum(150, norm, band_index, freqs, no_cycles, bands)