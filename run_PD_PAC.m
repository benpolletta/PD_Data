present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

for fo = 1:length(folders)
    
    folder = folders{fo};
   
    cd (folder)
    
    if ~exist('All_data.mat')
    
        PD_concatenate
    
    end
    
    if ~exist('All_data_dec.mat')
        
        PD_decimate
    
    else
        
        load('All_data_dec.mat','sampling_freq')
        
    end
       
    if ~exist('epoch1.txt')
        
        PD_epoch
    
    end
    
    if ~exist('PD_5min_master.list')
        
        PD_list
    
    end
        
    % challenge_list='PD_min_master.list';
    % 
    % challenge_descriptor='PD_min';
    % 
    % for m = -5:14
    %     labels{m+6} = num2str(m);
    % end
    
    challenge_list='PD_5min_master.list';
    
    challenge_descriptor='PD_5min';
    
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