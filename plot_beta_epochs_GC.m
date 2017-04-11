%run /project/crc-nak/brpp/startup

%cd /project/crc-nak/brpp/PD_Data/

subject_groups = {'st_m1_6OHDA', 'st_m1', 'st_stn'};
smooth_length = [5 20]*1000;
win_length = [333 1000];

for w = 1:2
    
    for s = 1:2
        
        for g = 1:3
            
            PD_beta_epochs_GC_collect([subject_groups{g}, '_subjects.mat'], 7, 2, win_length(w), smooth_length(s))
            
            PD_beta_epochs_GC_plot_group([subject_groups{g}, '_subjects.mat'], 7, 2, win_length(w), smooth_length(s))
            
        end
        
    end
    
end