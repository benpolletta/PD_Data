present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

for fo = 5:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
        
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
    
%     tic; wavelet_mouse_eeg_analysis_Jan_epsilon(sampling_freq,challenge_list,labels,{''},subplot_dims); toc;
    tic; wavelet_mouse_eeg_file_shuffle_IE(1000,0.99,challenge_list); toc;
    tic; wavelet_mouse_eeg_threshold_IE(1000,challenge_list,challenge_descriptor,labels,subplot_dims); toc;
    
    cd (present_dir)
    
end