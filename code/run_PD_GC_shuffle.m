function run_PD_GC_shuffle(no_shufs, epoch_length, maxorder)

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

sampling_freq = 1000;

f_length = sampling_freq*epoch_length + 1;

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
        
    listname = [prefix,'_channels_',num2str(epoch_length),'s_epochs.list'];
    
    epoch_list = text_read(listname,'%s%*[^\n]');
    
    no_epochs = length(epoch_list);
    
    All_GC_spec = nan(no_epochs,f_length,2);
    
    [striatal_indices,motor_indices] = random_pairs(no_shufs,no_epochs);
    
    no_shufs = length(striatal_indices);
    
    Orders = zeros(no_shufs,1);
    Errors = zeros(no_shufs,1);
    
    parfor s = 1:no_shufs
        
        local_epoch_list = epoch_list;
       
        striatal_name = local_epoch_list{striatal_indices(s)};
        motor_name = local_epoch_list{motor_indices(s)};
        
        s_data = load(striatal_name);
        s_data = s_data(:,1);
        
        m_data = load(motor_name);
        m_data = m_data(:,2);
        
        [moAIC,info,f] = mvgc_analysis([s_data m_data]',maxorder,'');
        
        Orders(s) = moAIC;
        
        Errors(s) = info.error;
        
        if info.error == 1
            
            display([striatal_name,', ',motor_name,' shuffle: ',info.errmsg])
            
        else
        
            GC_spec = reshape([reshape(f(2,1,:),1,f_length) reshape(f(1,2,:),1,f_length)],1,f_length,2);
            
            All_GC_spec(s,:,:) = GC_spec;
            
        end
        
    end
    
    save([listname(1:end-5),'_GCspec_shuffled.mat'],'Orders','Errors','All_GC_spec')
    
    cd (present_dir)
    
end