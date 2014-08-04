function run_PD_pre(outlier_check)

present_dir = pwd;

% load('initial_subjects.mat')

% load('st_m1_subjects.mat')

load('st_stn_subjects.mat')

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
%     PD_concatenate(prefix,1)
%     
%     PD_decimate(prefix,1)
    
    PD_concatenate_channels(prefix,1)
    
    PD_decimate_channels(prefix,1)
    
    if ~isempty(outlier_check)
        
        PD_epoch_list_artifacts(prefix,basetimes(fo),outlier_check)
        
        PD_plot_artifacts(prefix)
        
    end

    cd (present_dir)
    
end