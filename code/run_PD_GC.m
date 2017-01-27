function run_PD_GC(subjects_mat, peak_suffix, epoch_length, time_step, maxorder)

load(subjects_mat)

sampling_freq = 500;

present_dir = pwd;

f_length = sampling_freq*epoch_length + 1;

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    if isempty(time_step)
    
        listname = [prefix,'_channels_',num2str(epoch_length),'s_epochs.list'];
    
    else
       
        listname = [prefix,'_channels_',num2str(epoch_length),'s_by_',num2str(time_step),'s_epochs.list'];
        
    end
        
    epoch_list = text_read(listname,'%s%*[^\n]');
    
    no_epochs = length(epoch_list);
    
    All_GC_spec = nan(no_epochs,f_length,2);
    
    Orders = zeros(no_epochs,1);
    Errors = zeros(no_epochs,1);
    
    parfor e = 1:no_epochs
       
        epoch_name = epoch_list{e};
        
        if isempty(dir([epoch_name(1:end-4),'_GC.mat']))
            
            data = load(epoch_name);
            
            [moAIC,info,f] = mvgc_analysis(data',maxorder,epoch_name(1:end-4));
            
        else
            
            mvgc_results = load([epoch_name(1:end-4),'_GC.mat']);
            moAIC = mvgc_results.moAIC;
            info = mvgc_results.info;
            f = mvgc_results.f;
            
        end
        
        Orders(e) = moAIC;
        
        Errors(e) = info.error;
        
        if info.error == 1
            
            display([epoch_name,': ',info.errmsg])
            
        else
        
%             display(sprintf('%s: %d %d %d.',epoch_name,size(f)))
            
            GC_spec = reshape([reshape(f(2,1,:),1,f_length) reshape(f(1,2,:),1,f_length)],1,f_length,2);
            
            All_GC_spec(e,:,:) = GC_spec;
            
        end
        
    end
    
    save([listname(1:end-5),'_GCspec.mat'],'Orders','Errors','All_GC_spec')
    
    cd (present_dir)
    
end