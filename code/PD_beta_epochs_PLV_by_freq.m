function PD_beta_epochs_PLV_by_freq

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

f_bins = 7.5:5:32.5;

no_f_bins = length(f_bins);

freq_label = {'Striatum','Motor Ctx.'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    for ch = 1:2
         
        beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
        
        beta_P_listname = [beta_listname(1:end-5), 'P.list'];
       
        % Loading results of phase analyses.
        load([beta_P_listname(1:end-5),'byF.mat'])
        
        Good_MRV = nan(size(PbyF_mat, 1), no_f_bins - 1, 2);
        
        med_f = PbyF_mat(:, 1:(end-1), 4, :);        
        
        med_f = reshape(med_f, size(med_f, 1)*(no_f_bins - 1), 1, 2);
        
        med_f = reshape(med_f, size(med_f, 1), 2);
        
        for ch1 = 1:2
            
            for f = 1:(no_f_bins - 1)
                
                clear PLV_pval MRV n
                
                PLV_pval = PbyF_mat(:, f, 7, ch1);
                
                MRV = PbyF_mat(:, f, 5, ch1);
                
                n = PbyF_mat(:, f, 1, ch1);
                
                % Removing insignificant phase locking.
                Good_MRV_indicator = PLV_pval <= 0.01;
                
                Good_MRV(Good_MRV_indicator, f, ch1) = MRV(Good_MRV_indicator)./n(Good_MRV_indicator);
                
            end
          
        end
        
        Good_MRV = reshape(Good_MRV, size(Good_MRV, 1)*(no_f_bins - 1), 1, 2);
        
        Good_MRV = reshape(Good_MRV, size(Good_MRV, 1), 2);
                
        figure;
        
        for ch1 = 1:2
            
            %% Roseplot of PLV by freq.
            
            subplot(1,2,ch1)
            
            rose_plot(-angle(Good_MRV(:, ch1)), med_f(:, ch1), 20, f_bins);
            
            title({[folder,' ',freq_label{ch},' High Beta'];['Phase by ',freq_label{ch1},' \beta Freq.']})
            
        end
        
        save_as_pdf(gcf,[beta_listname(1:end-5),'_PLV_by_freq '])
        
    end
          
end