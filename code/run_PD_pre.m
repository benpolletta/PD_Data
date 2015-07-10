function run_PD_pre(subject_mat, format, detrend_window, outlier_check)

% detrend_window is in seconds.

present_dir = pwd;

load(subject_mat)

parfor fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    % PD_concatenate(prefix,1)
    % 
    % PD_decimate(prefix,1)
    
    PD_concatenate_channels(prefix, 1, format)
    
    PD_decimate_channels(prefix, detrend_window, 1)
    
    if ~isempty(outlier_check)
        
        PD_epoch_list_artifacts(prefix, basetimes(fo), outlier_check)
        
        PD_plot_artifacts(prefix)
        
    end

    cd (present_dir)
    
end
