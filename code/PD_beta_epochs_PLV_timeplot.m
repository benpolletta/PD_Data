function PD_beta_epochs_PLV_timeplot

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
        
        % Loading block numbers, epoch start and end indices.
        [blocks, epoch_starts, epoch_ends, ma_ch1, ma_ch2] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%f%f%*[^\n]');
        
        med_amp = [ma_ch1 ma_ch2];
        
        epoch_centers = round(mean([epoch_starts'; epoch_ends']));
        
        no_blocks = max(blocks);
        
        % Loading P_diff, F_smooth data.
        beta_pbf_name = [beta_listname(1:end-5),'_pbf_dp.txt'];
        
        beta_pbf_data = load(beta_pbf_name);
        
        blocks_by_dp = beta_pbf_data(:,1);
        
        All_Fs = beta_pbf_data(:,2:3);
        
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
        
        %% Plotting across blocks.
        
        figure;
        
        for ch1 = 1:2
            
            %% Roseplot of PLV by freq.
            
            subplot(2,2,ch1)
            
            rose_plot(-angle(Good_MRV), med_f(:, ch1), 20, 4:4:36);
            
            title({[folder,' ',freq_label{ch},' High Beta'];'Phase by \beta Freq.'})
            
            %% Roseplot of PLV by amp.
            
            subplot(2,2,2+ch1)
            
            rose_plot(-angle(Good_MRV), med_amp(:, ch1), 20, linspace(min(med_amp(:,ch1)), max(med_amp(:,ch1)), 10));
            
            title('Phase by \beta Amp.')
            
        end
        
        save_as_pdf(gcf,[beta_listname(1:end-5),'_PLV'])
        
        %% Plotting by block.
        
        for b = 1:no_blocks
            
            figure;
            
            block_indicator = blocks == b;
            
            block_start = min(epoch_starts(block_indicator));
            
            block_end = max(epoch_ends(block_indicator));
            
            t = (block_start:block_end)'/sampling_freq;
            
            Fs_block = All_Fs(blocks_by_dp == b,:);
            
            %% Plotting oscillations against GC.
            
            clear ax h2
            
            H_block = H(block_start:block_end,:,3);
            
            block_epochs = block_start < epoch_centers & epoch_centers < block_end;
            
            subplot(4,1,1)
            
            [ax, ~, h2] = plotyy(t,real(H_block),epoch_centers(block_epochs)/sampling_freq, -angle(Good_MRV(block_epochs)));
            
            title('\beta Oscillation and Phase-Locking')
            
            xlabel('Time (s)')
            
            axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
            
            set(ax,{'ycolor'},{'black';'black'})
            
            set(get(ax(1),'YLabel'),'String','\beta Oscillation')
            
            legend(ax(1), freq_label,'Location','NorthWest')
            
            set(get(ax(2),'YLabel'),'String','Significant Phase-Locking')
            
            set(h2,'LineStyle','*')
            
            %% Plotting freq. against GC.
            
            clear ax h2
            
            subplot(4,1,2)
            
            [ax, ~, h2] = plotyy(t,Fs_block,epoch_centers(block_epochs)/sampling_freq, -angle(Good_MRV(block_epochs)));
            
            title('\beta Frequency and Phase Locking')
            
            xlabel('Time (s)')
            
            axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
            
            set(ax,{'ycolor'},{'black';'black'})
            
            set(get(ax(1),'YLabel'),'String','\beta Frequency')
            
            set(get(ax(2),'YLabel'),'String','Significant Phase-Locking')
            
            set(h2,'LineStyle','*')
            
            %% Plotting amp. against GC.
            
            clear ax h2
            
            A_block = A(block_start:block_end,:,3);
            
            subplot(4,1,3)
            
            [ax, ~, h2] = plotyy(t,A_block,epoch_centers(block_epochs)/sampling_freq, -angle(Good_MRV(block_epochs)));
            
            title('\beta Amplitude and Phase Locking')
            
            xlabel('Time (s)')
            
            axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
            
            set(ax,{'ycolor'},{'black';'black'})
            
            set(get(ax(1),'YLabel'),'String','\beta Amplitude')
            
            set(get(ax(2),'YLabel'),'String','Significant Phase-Locking')
            
            set(h2,'LineStyle','*')
            
            for ch1 = 1:2
                
                %% Roseplot of PLV by freq.
                
                subplot(4, 4, 12 + ch1)
                
                rose_plot(-angle(Good_MRV(block_epochs)), med_f(block_epochs, ch1), 20, 4:4:36);
                
                title({[folder,' ',freq_label{ch},' High Beta Block ',num2str(b)];'Phase by \beta Freq.'})
                
                %% Roseplot of PLV by amp.
                
                subplot(4, 4, 14 + ch1)
                
                rose_plot(-angle(Good_MRV(block_epochs)), med_amp(block_epochs, ch1), 20, linspace(min(med_amp(block_epochs,ch1)), max(med_amp(block_epochs,ch1)), 10));
                
                title({[folder,' ',freq_label{ch},' High Beta Block ',num2str(b)];'Phase by \beta Amp.'})
                
            end
            
            save_as_pdf(gcf,[beta_listname(1:end-5),'_block',num2str(b),'_PLV'])
            
            close(gcf)
            
        end
        
    end
          
end