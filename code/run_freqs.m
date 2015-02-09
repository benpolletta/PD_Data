function run_freqs

dirs = {'.', '.', '.', 'Carb_Opto', '6OHDA_Opto'};

subject_matnames = {'st_m1', 'st_stn', 'st_m1_6OHDA', 'st_m1_emxarch', 'st_m1_6OHDA_opto'};

norms = {'', '_pct', '_zscore'};

freqs = {[], [8:.5:30], [8:.5:30]};

no_cycles = {[], linspace(3, 7, length(freqs{2})), linspace(7, 7^2, length(freqs{3}))};

bands = {[], [8 13; 13 18; 18 25; 13 25; 18 30], [8 13; 13 18; 18 25; 13 25; 18 30]};

for k = 1:length(freqs)
    
    for i = 1:length(subject_matnames)
        
        cd (dirs{i})
        
        for j = 1:length(norms)
            
            beta_blocks_rel_infusion_freqs([subject_matnames{i}, '_subjects.mat'], norms{j}, freqs{k}, no_cycles{k}, bands{k})
            
        end
        
        cd '/home/bp/PD_Data'
        
    end
    
    for i = 1:length(subject_matnames)
        
        cd (dirs{i})
        
        for j = 1:length(norms)
            
            beta_blocks_rel_infusion_freq_plot([subject_matnames{i}, '_subjects.mat'], norms{j}, '', freqs{k}, no_cycles{k}, bands{k})
            
        end
        
        cd '/home/bp/PD_Data'
        
    end
    
end