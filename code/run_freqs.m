function run_freqs

home_dir = pwd;

dirs = {'.', '.', '.', 'Carb_Opto', '6OHDA_Opto'};

subject_matnames = {'st_m1', 'st_stn', 'st_m1_6OHDA'}; %, 'st_m1_emxarch', 'st_m1_6OHDA_opto'};

norms = {'', '_pct', '_zscore'};

band_indices = [3 6 6];

freqs = {[], [8:.5:30], [8:.5:30]};

no_cycles = {[], linspace(3, 7, length(freqs{2})), linspace(7, 7^2, length(freqs{3}))};

bands = {[], [8 13; 13 18; 18 25; 13 25; 18 30; 8 30], [8 13; 13 18; 18 25; 13 25; 18 30; 8 30]};

epoch_secs = {150, []};

for k = 1:length(freqs)
    
    for i = 1:length(subject_matnames)
    
        cd (dirs{i})
    
        for j = 1:length(norms)
            
            for l = 1:length(epoch_secs)
                
                beta_blocks_rel_infusion_freqs([subject_matnames{i}, '_subjects.mat'], norms{j}, band_indices(k), epoch_secs{l}, freqs{k}, no_cycles{k}, bands{k})
                
                beta_blocks_rel_infusion_freq_plot([subject_matnames{i}, '_subjects.mat'], norms{j}, '', band_indices(k), epoch_secs{l}, freqs{k}, no_cycles{k}, bands{k})
                
            end
            
        end
        
        cd (home_dir)
        
    end
    
end
