present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    if isempty(dir([prefix,'_all_data.mat']))
    
        PD_concatenate(prefix)
    
    end
    
    if isempty(dir([prefix,'_all_data_dec.mat']))
        
        PD_decimate(prefix)
        
    end
        
    load([prefix,'_all_data_dec.mat'],'sampling_freq')
       
    if isempty(dir([prefix,'_epoch1.txt']))
        
        PD_epoch(prefix)
    
    end
    
    if isempty(dir([prefix,'_5min_master.list']))
        
        PD_list(prefix)
    
    end
        
    % challenge_list=[prefix,'_min_master.list'];
    % 
    % challenge_descriptor=[prefix,'_min'];
    % 
    % for m = -5:14
    %     labels{m+6} = num2str(m);
    % end
    
    challenge_list=[prefix,'_5min_master.list'];
    
    peak_averaged_signal_batch_condition_parallel(challenge_list,140,4,4,sampling_freq);
    
    cd (present_dir)
    
end