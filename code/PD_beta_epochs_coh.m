function PD_beta_epochs_coh(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

pd_label = {'pre', 'post'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    for ch = 1:no_channels
        
        for pd = 1:2
            
            beta_listname = [subj_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,'.list'];
            % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
            
            beta_list = text_read(beta_listname, '%s%*[^\n]');
            
            no_epochs = length(beta_list);
            
            All_coh = nan(no_epochs, window_length);
            
            for e = 1:no_epochs
                
                data_name = beta_list{e};
                
                data = load(data_name);
                
                data_hat = fft(data);
                
                data_xspec = data_hat(:, 1) .* conj(data_hat(:, 2));
                
                data_xspec_norm = data_xspec * diag(1 ./ sum(data_hat .* conj(data_hat)));
                
                All_coh(e, :) = data_xspec_norm';
                
            end
            
            save([beta_listname(1:end-5),'_coh.mat'], 'All_coh')
            
        end
        
    end
    
end