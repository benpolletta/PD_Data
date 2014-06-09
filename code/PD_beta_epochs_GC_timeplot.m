function PD_beta_epochs_GC_timeplot

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

f_bins = 7.5:5:32.5;

no_f_bins = length(f_bins);

sampling_freq = 1000;
            
freq_label = {'Striatum','Motor Ctx.'};

gc_label = {'Str. --> Motor','Motor --> Str.'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
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
       
        % Loading results of GC and phase analyses.
        load([beta_listname(1:end-5),'GC.mat'])
        
        load([beta_P_listname(1:end-5),'byF.mat'])
        
        Good_GC = nan(size(All_GC,1),2);
        
        for ch1 = 1:2
            
            % Removing insignificant or erroneous GC.
            GC = All_GC(:, 2 + ch1);
            
            GC_pval = All_GC(:, 4 + ch1);
            
            Good_GC_indicator = GC_pval <= 0.01 & ~isnan(GC);
            
            Good_GC(Good_GC_indicator, ch1) = GC(Good_GC_indicator);
            
        end
        
        %% Plotting across blocks.
        
        figure;
        
        for ch1 = 1:2
            
            %% Correlation of GC by freq.
            
            subplot(2,4,ch1)
            
            r = nancorr(PbyF_mat(:, no_f_bins, 4, ch1), Good_GC(:, ch1));
            
            [hist, centers] = hist3([PbyF_mat(:, no_f_bins, 4, ch1) Good_GC(:, ch1)], [25 25]);
            
            colormap('default')
            
            imagesc(centers{1},centers{2},hist/sum(sum(hist)))
            
            h = colorbar;
            
            title({[folder,' ',freq_label{ch},' High Beta'];['Correlation = ',num2str(r)]})
            
            ylabel(h,'Prop. Observed')
            
            axis xy
            
            xlabel(['Median Freq. (',freq_label{ch1},')'])
            
            ylabel(['GC (',gc_label{ch1},')'])
            
            subplot(2,4,2+ch1)
            
            r = nancorr(PbyF_mat(:, no_f_bins, 4, 3-ch1), Good_GC(:, ch1));
            
            [hist, centers] = hist3([PbyF_mat(:, no_f_bins, 4, 3-ch1) Good_GC(:, ch1)], [25 25]);
            
            colormap('default')
            
            imagesc(centers{1},centers{2},hist/sum(sum(hist)))
            
            h = colorbar;
            
            title({[folder,' ',freq_label{ch},' High Beta'];['Correlation = ',num2str(r)]})
            
            ylabel(h,'Prop. Observed')
            
            axis xy
            
            xlabel(['Median Freq. (',freq_label{3-ch1},')'])
            
            ylabel(['GC (',gc_label{ch1},')'])
            
            % Correlation of GC by amp.
            
            subplot(2,4,4+ch1)
            
            r = nancorr(med_amp(:, ch1), Good_GC(:, ch1));
            
            [hist, centers] = hist3([med_amp(:, ch1) Good_GC(:, ch1)], [25 25]);
            
            colormap('default')
            
            imagesc(centers{1},centers{2},hist/sum(sum(hist)))
            
            title(['Correlation = ',num2str(r)])
            
            h = colorbar;
            
            ylabel(h,'Prop. Observed')
            
            axis xy
            
            xlabel(['Median Amp. (',freq_label{ch1},' \beta)'])
            
            ylabel(['GC (',gc_label{ch1},')'])
            
            subplot(2,4,6+ch1)
            
            r = nancorr(med_amp(:, 3-ch1), Good_GC(:, ch1));
            
            [hist, centers] = hist3([med_amp(:, 3-ch1) Good_GC(:, ch1)], [25 25]);
            
            colormap('default')
            
            imagesc(centers{1},centers{2},hist/sum(sum(hist)))
            
            title(['Correlation = ',num2str(r)])
            
            h = colorbar;
            
            ylabel(h,'Prop. Observed')
            
            axis xy
            
            xlabel(['Median Amp. (',freq_label{3-ch1},' \beta)'])
            
            ylabel(['GC (',gc_label{ch1},')'])
            
        end
        
        save_as_pdf(gcf,[beta_listname(1:end-5),'_GC'])
        
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
            
            subplot(5,1,1)
            
            [ax, ~, h2] = plotyy(t,real(H_block),epoch_centers(block_epochs)/sampling_freq,Good_GC(block_epochs,:));
            
            title('\beta Oscillation and Granger Causality')
            
            xlabel('Time (s)')
            
            axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
            
            set(ax,{'ycolor'},{'black';'black'})
            
            set(get(ax(1),'YLabel'),'String','\beta Oscillation')
            
            legend(ax(1), freq_label,'Location','NorthWest')
            
            set(get(ax(2),'YLabel'),'String','Significant GC')
            
            legend(ax(2), gc_label,'Location','NorthEast')
            
            set(h2,'LineStyle','*')
            
            %% Plotting freq. against GC.
            
            clear ax h2
            
            subplot(5,1,2)
            
            [ax, ~, h2] = plotyy(t,Fs_block,epoch_centers(block_epochs)/sampling_freq,Good_GC(block_epochs,:));
            
            title('\beta Frequency and Granger Causality')
            
            xlabel('Time (s)')
            
            axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
            
            set(ax,{'ycolor'},{'black';'black'})
            
            set(get(ax(1),'YLabel'),'String','\beta Frequency')
            
            set(get(ax(2),'YLabel'),'String','Significant GC')
            
            set(h2,'LineStyle','*')
            
            %% Plotting amp. against GC.
            
            clear ax h2
            
            A_block = A(block_start:block_end,:,3);
            
            subplot(5,1,3)
            
            [ax, ~, h2] = plotyy(t,A_block,epoch_centers(block_epochs)/sampling_freq,Good_GC(block_epochs,:));
            
            title('\beta Amplitude and Granger Causality')
            
            xlabel('Time (s)')
            
            axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
            
            set(ax,{'ycolor'},{'black';'black'})
            
            set(get(ax(1),'YLabel'),'String','\beta Amplitude')
            
            set(get(ax(2),'YLabel'),'String','Significant GC')
            
            set(h2,'LineStyle','*')
                
            for ch1 = 1:2
            
                %% Correlation of GC by freq.
                
                subplot(5,4,12+ch1)
                
                [hist, centers] = hist3([PbyF_mat(block_epochs, no_f_bins, 4, ch1) Good_GC(block_epochs, ch1)],[20 20]);
                
                colormap('default')
                
                imagesc(centers{1},centers{2},hist/sum(sum(hist)))
                
                h = colorbar;
                
                ylabel(h,'Prop. Observed')
                
                axis xy
                
                xlabel(['Median Freq. (',freq_label{ch1},')'])
                
                ylabel(['GC (',gc_label{ch1},')'])
                
                subplot(5,4,14+ch1)
                
                [hist, centers] = hist3([PbyF_mat(block_epochs, no_f_bins, 4, 3-ch1) Good_GC(block_epochs, ch1)],[20 20]);
                
                colormap('default')
                
                imagesc(centers{1},centers{2},hist/sum(sum(hist)))
                
                h = colorbar;
                
                ylabel(h,'Prop. Observed')
                
                axis xy
                
                xlabel(['Median Freq. (',freq_label{3-ch1},')'])
                
                ylabel(['GC (',gc_label{ch1},')'])
            
                %% Correlation of GC by amp.
                
                subplot(5,4,16+ch1)
                
                [hist, centers] = hist3([med_amp(block_epochs, ch1) Good_GC(block_epochs, ch1)], [20 20]);
                
                colormap('default')
                
                imagesc(centers{1},centers{2},hist/sum(sum(hist)))
                
                h = colorbar;
                
                ylabel(h,'Prop. Observed')
                
                axis xy
                
                xlabel(['Median Amp. (',freq_label{ch1},' \beta)'])
                
                ylabel(['GC (',gc_label{ch1},')'])
                
                subplot(5,4,18+ch1)
                
                [hist, centers] = hist3([med_amp(block_epochs, 3-ch1) Good_GC(block_epochs, ch1)], [20 20]);
                
                colormap('default')
                
                imagesc(centers{1},centers{2},hist/sum(sum(hist)))
                
                h = colorbar;
                
                ylabel(h,'Prop. Observed')
                
                axis xy
                
                xlabel(['Median Amp. (',freq_label{3-ch1},' \beta)'])
                
                ylabel(['GC (',gc_label{ch1},')'])
                
            end
            
            save_as_pdf(gcf,[beta_listname(1:end-5),'_block',num2str(b),'_GC'])
            
            close(gcf)
            
        end
        
    end
          
end