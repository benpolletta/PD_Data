function PD_beta_epochs_PbyF

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

varnames = {'f_datapoints','f_mean','f_std','mrv','ph_zstat','ph_pval'};

f_bins = 7.5:5:32.5;
            
no_f_bins = length(f_bins)-1; % Finding number of frequency bins.

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    for ch = 1:2
         
        beta_listname = [subj_name,'_ch',num2str(ch), '_betaP.list'];
        
        beta_list = text_read(beta_listname, '%s%*[^\n]');
        
        no_epochs = length(beta_list);
        
        PbyF_mat = nan(no_epochs, no_f_bins, 6, 2);
        
        parfor e = 1:no_epochs
           
            P_name = beta_list{e};
            
            pbf = PbyF_epoch([P_name(1:end-5),'P.txt'], f_bins);

            PbyF_mat(e, :, :, :) = pbf;
            
        end
        
        save([beta_listname(1:end-5),'byF.mat'],'PbyF_mat','varnames')
        
    end
    
end