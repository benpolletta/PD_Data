function run_PD_GC_shuffle(no_shufs)

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

sampling_freq = 1000;

epoch_length = 5*sampling_freq;

f_length = epoch_length/2 + 1;

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
        
    epoch_list = text_read([prefix,'_channels_epochs.list'],'%s%*[^\n]');
    
    no_epochs = length(epoch_list);
    
    All_GC_spec = nan(no_epochs,f_length,2);
    
    [striatal_indices,motor_indices] = random_pairs(no_shufs,no_epochs);
    
    no_shufs = length(striatal_indices);
    
    parfor s = 1:no_shufs
        
        local_epoch_list = epoch_list;
       
        striatal_name = local_epoch_list{striatal_indices(s)};
        motor_name = local_epoch_list{motor_indices(s)};
        
        s_data = load(striatal_name);
        s_data = s_data(:,1);
        
        m_data = load(motor_name);
        m_data = m_data(:,2);
        
        [moAIC,f] = mvgc_analysis([s_data m_data]',100,'');
        
        Orders(s) = moAIC;
        
        GC_spec = reshape([reshape(f(2,1,:),1,f_length) reshape(f(1,2,:),1,f_length)],1,f_length,2);
        
        All_GC_spec(s,:,:) = GC_spec;
        
    end
    
    save([prefix,'_channels_epochs_GCspec_shuffled.mat'],'Orders','All_GC_spec')
    
    cd (present_dir)
    
end