function run_PD_PbyF(subject_mat)

present_dir = pwd;

load(subject_mat)

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    PD_concatenate_channels(prefix,1)
    
    PD_decimate_channels(prefix,1)

    cd (present_dir)
    
end

PD_bandpass_channels(subject_mat)

PD_beta_epochs(subject_mat,2,round(1000/3),100)

PD_beta_epochs_roseplot_by_datapoint(subject_mat)

PD_beta_epochs_rel_infusion(subject_mat,2,round(1000/3))

PD_beta_epochs_rel_infusion_roseplot_by_datapoint(subject_mat)