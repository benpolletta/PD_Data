function run_fix

subject_matname = {'st_stn', 'st_m1'};

freqs = {[], [8:.5:30], [8:.5:30]};

no_cycles = {[], linspace(3, 7, length(freqs{2})), linspace(7, 7^2, length(freqs{3}))};

bands = {[], [8 13; 13 18; 18 25; 13 25; 18 30], [8 13; 13 18; 18 25; 13 25; 18 30]};

new_bands = {[], [8 13; 13 18; 18 25; 13 25; 18 30; 8 30], [8 13; 13 18; 18 25; 13 25; 18 30; 8 30]};

for i = 1:2
    
    for j = 1:3
        
        if j == 1 && i == 2
            
            PD_fix_Spec_pct([subject_matname{i}, '_subjects.mat'], freqs{j}, no_cycles{j}, bands{j})
            
        end
        
        PD_fix_wav_BP([subject_matname{i}, '_subjects.mat'], freqs{j}, no_cycles{j}, bands{j}, new_bands{j})
        
    end
    
end