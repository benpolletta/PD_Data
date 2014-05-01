function run_PD_GC

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

sampling_freq = 1000;

epoch_length = 5*1000;

f_length = epoch_length/2 + 1;

f = (sampling_freq/2)*(1:f_length)/f_length;

for fo = 5:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
        
    epoch_list = text_read([prefix,'_channels_epochs.list'],'%s%*[^\n]');
    
    no_epochs = length(epoch_list);
    
    All_GC_spec = nan(f_length,no_epochs,2);
    
    for e = 1:no_epochs
       
        epoch_name = epoch_list{e};
        
        data = load(epoch_name);
        
        [moAIC,f] = mvgc_analysis(data',100,epoch_name(1:end-4));
        
        Orders(e) = moAIC;
        
        All_GC_spec(:,e,1) = reshape(f(2,1,:),f_length,1);
        All_GC_spec(:,e,2) = reshape(f(1,2,:),f_length,1);
        
    end
    
    save([prefix,'_channels_epochs_all_GCspec.mat'],'Orders','All_GC_spec')
    
    cd (present_dir)
    
end