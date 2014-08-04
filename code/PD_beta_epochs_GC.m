function PD_beta_epochs_GC(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_index = {1, 2, 1:2, 1:2, 1, 2, 1, 2};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

pd_label = {'pre', 'post'};

list_suffix = {'', '_A'};

no_lists = length(list_suffix);

f_length = win_size + 2;

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    for ch = 1:no_channels
        
        for pd = 1:2
            
            for list = 1:no_lists
                
                beta_listname = [subj_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,list_suffix{list},'.list'];
                % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
                
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
                
                save([beta_listname(1:end-5),'_GC.mat'],'All_GC','All_GC_spec')
                
            end
            
        end
        
    end
    
end