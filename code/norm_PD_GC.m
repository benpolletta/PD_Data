function norm_PD_GC(epoch_length)

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    %% Shuffles.
    
    shuf_name = [prefix,'_channels_',num2str(epoch_length),'s_epochs_GCspec_shuffled.mat'];
    
    load(shuf_name)
    
    All_GC_spec = abs(All_GC_spec);
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
    shuffled_mean = nanmean(All_GC_spec);
    
    shuffled_std = nanstd(All_GC_spec);
    
    save([shuf_name(1:end-4),'_stats.mat'],'shuffled_mean','shuffled_std')
    
    %% Observed.
    
    clear All_GC_spec
    
    obs_name = [prefix,'_channels_',num2str(epoch_length),'s_epochs_GCspec.mat'];
    
    load(obs_name)
    
    All_GC_spec = abs(All_GC_spec);
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
    spec_mean = nanmean(All_GC_spec);
    
    spec_std = nanstd(All_GC_spec);
    
    % Computing normalized GC.
    
    All_GC_norm = nan(size(All_GC_spec));
    
    for i=1:3
        
        All_GC_norm(:,:,i) = (All_GC_spec(:,:,i) - ones(size(All_GC_spec(:,:,i)))*diag(spec_mean(:,:,i)))./(ones(size(All_GC_spec(:,:,i)))*diag(spec_std(:,:,i)));
        
    end
    
    All_GC_zs = nan(size(All_GC_spec));
    
    for i=1:3
        
        All_GC_zs(:,:,i) = (All_GC_spec(:,:,i) - ones(size(All_GC_spec(:,:,i)))*diag(shuffled_mean(:,:,i)))./(ones(size(All_GC_spec(:,:,i)))*diag(shuffled_std(:,:,i)));
        
    end
    
    All_GC_zs(:,:,4) = -diff(All_GC_zs(:,:,1:2),[],3);
    
    save([obs_name(1:end-4),'_thresh.mat'],'All_GC_norm','All_GC_zs')
    
    zs_mean = nanmean(All_GC_zs);
    
    zs_std = nanstd(All_GC_zs);
    
    save([obs_name(1:end-4),'_stats.mat'],'spec_mean','spec_std','zs_mean','zs_std')
    
    cd (present_dir)
    
end