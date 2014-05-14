function run_PD_GC_by_phase

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    cd (folder)
    
    master_listname = [folder,'_all_channel_data_dec_phase_master.list'];
    
    master_list = text_read(master_listname,'%s%*[^\n]');
    
    no_lists = length(master_list);
    
    for l = 1:no_lists
        
        listname = master_list{l};
        
        epoch_list = text_read(listname,'%s%*[^\n]');
        
        no_epochs = length(epoch_list);
        
        All_GC = nan(no_epochs,2);
        All_pval = nan(no_epochs,2);
        All_sig = nan(no_epochs,2);
        
        Orders = zeros(no_epochs,1);
        Errors = zeros(no_epochs,1);
        
        parfor e = 58:no_epochs
            
            epoch_name = epoch_list{e};
            
            epoch_name = epoch_name(8:end);
            
            if isempty(dir([epoch_name(1:end-4),'_GC.mat']))
                
                data = load(epoch_name);
                
                data = reshape(data',2,length(data)/2);
                
                [moAIC,info,F,pval,sig] = mvgc_analysis(data',20,epoch_name(1:end-4),0);
                
            else
                
                mvgc_results = load([epoch_name(1:end-4),'_GC.mat']);
                moAIC = mvgc_results.moAIC;
                info = mvgc_results.info;
                F = mvgc_results.F;
                pval = mvgc_results.pval;
                sig = mvgc_results.sig;
                
            end
            
            Orders(e) = moAIC;
            
            Errors(e) = info.error;
            
            if info.error == 1
                
                display([epoch_name,': ',info.errmsg])
                
            else
                
                All_GC(e,:) = diag(flipud(F));
                
                All_pval(e,:) = diag(flipud(pval));
                
                All_sig(e,:) = diag(flipud(sig));
                
            end
            
        end
        
        save([listname(1:end-5),'_GC.mat'],'Orders','Errors','All_GC','All_pval','All_sig')
        
        cd (present_dir)
        
    end
    
end