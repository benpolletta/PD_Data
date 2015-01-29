function run_freqs

dirs = {'.', '.', '.', 'Carb_Opto', '6OHDA_Opto'};

subject_matnames = {'st_m1', 'st_stn', 'st_m1_6OHDA', 'st_m1_emxarch', 'st_m1_6OHDA_opto'};

norms = {'', '_pct', '_zscore'};

for i = 1:length(subject_matnames)
    
    cd (dirs{i})
    
    for j = 2 % 1:length(norms)
        
        beta_blocks_rel_infusion_freqs([subject_matnames{i}, '_subjects.mat'], norms{j})
        
        % beta_blocks_rel_infusion_freq_plot([subject_matnames{i}, '_subjects.mat'], norms{j}, '')
        
    end
    
    cd '/home/bp/PD_Data'
    
end

for i = 1:length(subject_matnames)
    
    cd (dirs{i})
    
    for j = 1:length(norms)
        
        % beta_blocks_rel_infusion_freqs([subject_matnames{i}, '_subjects.mat'], norms{j})
        
        beta_blocks_rel_infusion_freq_plot([subject_matnames{i}, '_subjects.mat'], norms{j}, '')
        
    end
    
    cd '/home/bp/PD_Data'
    
end