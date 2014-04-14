present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

for fo = 1:1%:length(folders)
    
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
    
    challenge_descriptor=[prefix,'_5min'];
    
    for m = -1:2
        labels{m+2} = ['5 Min. Period Rel. Injection: ',num2str(m),'.'];
    end
    
    [sr, sc] = subplot_size(length(labels));
    subplot_dims = [sr, sc];
    
    tic; wavelet_mouse_eeg_analysis_Jan_epsilon(sampling_freq,challenge_list,labels,{''},subplot_dims); toc;
    tic; wavelet_mouse_eeg_file_shuffle_IE(1000,0.99,challenge_list); toc;
    tic; wavelet_mouse_eeg_threshold_IE(1000,challenge_list,challenge_descriptor,labels,subplot_dims); toc;
    
    cd (present_dir)
    
end