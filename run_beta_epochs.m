% run ../startup
% 
% cd /project/crc-nak/brpp/PD_Data/

subject_names = {'st_m1', 'st_stn', 'st_m1_60HDA'};

win_pairs = [333 5000; 1000 5000; 333 20000; 2000 20000];

for s = 1:3
    
    for w = 1:4
   
        PD_beta_epochs_rel_infusion([subject_names{s},'_subjects.mat'], 7, 2, win_pairs(w, 1), win_pairs(w, 2));
        
    end
    
end