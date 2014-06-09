function PD_beta_epochs_PLV_by_GC

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

f_bins = 7.5:5:32.5;

no_f_bins = length(f_bins);
            
freq_label = {'Striatum','Motor Ctx.'};
            
gc_label = {'Str. --> Motor','Motor --> Str.'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    for ch = 1:2
         
        beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
        
        beta_P_listname = [beta_listname(1:end-5), 'P.list'];
       
        % Loading results of phase analyses.
        load([beta_P_listname(1:end-5),'byF.mat'])
        
        clear PLV_pval MRV n med_f
        
        PLV_pval = PbyF_mat(:, no_f_bins, 7, 1);
        
        MRV = PbyF_mat(:, no_f_bins, 5, 1);
        
        n = PbyF_mat(:, no_f_bins, 1, 1);
        
        med_f(:,1) = PbyF_mat(:, no_f_bins, 4, 1); med_f(:,2) = PbyF_mat(:, no_f_bins, 4, 2);
        
        % Removing insignificant phase locking.
        Good_MRV_indicator = PLV_pval <= 0.01;
        
        Good_MRV = nan(size(MRV));
        
        Good_MRV(Good_MRV_indicator) = MRV(Good_MRV_indicator)./n(Good_MRV_indicator);
       
        % Loading results of GC analysis.
        load([beta_listname(1:end-5),'GC.mat'])
        
        Good_GC = nan(size(All_GC,1),2);
            
        % Removing insignificant or erroneous GC.
        for ch1 = 1:2
            
            GC = All_GC(:, 2 + ch1);
            
            GC_pval = All_GC(:, 4 + ch1);
            
            Good_GC_indicator = GC_pval <= 0.01 & ~isnan(GC);
            
            Good_GC(Good_GC_indicator, ch1) = GC(Good_GC_indicator);
            
        end
        
        %% Plotting across blocks.
        
        figure;
        
        for ch1 = 1:2
            
            %% Roseplot of PLV by GC.
            
            subplot(1,2,ch1)
            
            rose_plot(-angle(Good_MRV), Good_GC(:, ch1), 20, linspace(min(Good_GC(:,ch1)), max(Good_GC(:,ch1)), 10));
            
            title({[folder,' High ',freq_label{ch},' Beta'];['Phase by GC (',gc_label{ch1},').']})
            
        end
        
        save_as_pdf(gcf,[beta_listname(1:end-5),'_PLV_by_GC'])
        
    end
          
end