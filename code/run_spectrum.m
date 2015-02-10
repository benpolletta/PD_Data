function run_spectrum

subject_matname = {'st_stn', 'st_m1', 'st_m1_6OHDA'};

band_indices = [3 6 6];

freqs = {[], [8:.5:30], [8:.5:30]};

no_cycles = {[], linspace(3, 7, length(freqs{2})), linspace(7, 7^2, length(freqs{3}))};

bands = {[], [8 13; 13 18; 18 25; 13 25; 18 30], [8 13; 13 18; 18 25; 13 25; 18 30]};

for i = 1 % 1:2
    
    for j = 1:3
        
        PD_beta_blocks_rel_infusion_pre_post_spectrum([subject_matname{i}, '_subjects.mat'], 150, '_by_STR', '', band_indices(j), freqs{j}, no_cycles{j}, bands{j})
        
        PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual([subject_matname{i}, '_subjects.mat'], 150, '_by_STR', '', band_indices(j), freqs{j}, no_cycles{j}, bands{j})
        
    end
    
end