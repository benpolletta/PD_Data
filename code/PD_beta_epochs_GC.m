function PD_beta_epochs_GC(win_length)

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

f_length = win_length + 2;

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    for ch = 1:2
         
        beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
        
        beta_list = text_read(beta_listname, '%s%*[^\n]');
        
        no_epochs = length(beta_list);
        
        All_GC = nan(no_epochs, 8);
        
        All_GC_spec = nan(no_epochs, f_length, 2);
        
        parfor e = 1:no_epochs
           
            data_name = beta_list{e};
            
            data = load(data_name);
            
            [moAIC, info, F, pval, sig, f] = mvgc_analysis(data, 20, data_name(1:end-4), 2);
            
            All_GC(e, :) = [moAIC info.error diag(flipud(F))' diag(flipud(pval))' diag(flipud(sig))'];
            
            All_GC_spec(e, :, :) = reshape([reshape(f(2,1,:),1,f_length) reshape(f(1,2,:),1,f_length)],1,f_length,2);
            
        end
        
        save([beta_listname(1:end-5),'GC.mat'],'All_GC','All_GC_spec')
        
    end
    
end