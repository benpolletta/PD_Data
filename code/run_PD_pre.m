function run_PD_pre(outlier_check)

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

basetimes = [300 1200 1800 600 1800];

infusetimes = [390 240 300 450 510];

for fo = 5:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
%     PD_concatenate(prefix,1)
    
%     PD_decimate(prefix,1)
    
    PD_epoch_list_artifacts(prefix,basetimes(fo),outlier_check)
    
    PD_plot_artifacts(prefix)

    cd (present_dir)
    
end