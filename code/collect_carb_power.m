function collect_carb_power(norm, freqs, no_cycles, bands)


matnames = {'st_m1', 'st_stn', 'st_m1_6OHDA'};

period_labels = {'_by_STR', '_by_STR', ''};

matlabpool open 3

parfor i = 1:3

    PD_beta_blocks_rel_infusion_pre_post_power([matnames{i}, '_subjects.mat'],...
        150, period_labels{i}, norm, freqs, no_cycles, bands)
    
end

collect_striatal_power_density(150, '_power', norm, freqs, no_cycles, bands)

combine_carb_6OHDA_power_density(150, '_power', norm, freqs, no_cycles, bands)

matlabpool close
