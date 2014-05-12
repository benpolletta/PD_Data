function run_GC_pre(epoch_length,time_step,outlier_check)

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

basetimes = [300 1200 1800 600 1800];

infusetimes = [390 240 300 450 510];

load('BetaTimes')

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
%     PD_concatenate_channels(prefix,1)
    
%     PD_decimate_channels(prefix,1)

    if isempty(time_step)

        PD_epoch_list_artifacts_channels(prefix, epoch_length, outlier_check)

    else

        PD_epoch_list_artifacts_channels_sliding(prefix, epoch_length, time_step, outlier_check)

    end
    
    cd (present_dir)
    
end

% get_beta