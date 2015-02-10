function run_fix

home_dir = pwd;

dir = {'', '', '', 'Carb_Opto', '6OHDA_Opto'};

subject_matname = {'st_stn', 'st_m1', 'st_m1_6OHDA', 'st_m1_emxarch', 'st_m1_6OHDA_opto'};

freqs = {[], [8:.5:30], [8:.5:30]};

no_cycles = {[], linspace(3, 7, length(freqs{2})), linspace(7, 7^2, length(freqs{3}))};

bands = {[], [8 13; 13 18; 18 25; 13 25; 18 30], [8 13; 13 18; 18 25; 13 25; 18 30]};

new_bands = {[], [8 13; 13 18; 18 25; 13 25; 18 30; 8 30], [8 13; 13 18; 18 25; 13 25; 18 30; 8 30]};

for i = 4:5 % 1:5
    
    for j = 2:3
        
        % PD_fix_Spec_pct([subject_matname{i}, '_subjects.mat'], freqs{j}, no_cycles{j}, bands{j})
        
        if i == 1
            
            if j > 1
                
                PD_fix_wav_BP([subject_matname{i}, '_subjects.mat'], freqs{j}, no_cycles{j}, bands{j}, new_bands{j})
                
            end
            
        else
            
            cd (dir{i})
            
            PD_fix_wav_BP([subject_matname{i}, '_subjects.mat'], freqs{j}, no_cycles{j}, bands{j}, new_bands{j})
            
            cd (home_dir)
            
        end
            
    end
    
end