function beta_blocks_consolidated_plot(subject_mat, time_window, percent)

PD_struct = PD_initialize(subject_mat);

for fo = 1:PD_struct.no_folders
    
    folder = PD_struct.folders{fo};
    
    prefix = PD_struct.prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    cons_name = [subj_name, '_2sd_BP_high_', num2str(time_window/PD_struct.sampling_freq), 's_', num2str(percent), 'pct_consolidated'];
    
    cons_title = [folder, ', High Power, ', num2str(percent), ' Percent of ', num2str(time_window/PD_struct.sampling_freq), ' s'];
    
    load([cons_name, '.mat'], 'BP_high_consolidated')
    
    load([subj_name, '_all_channel_data_dec.mat'])
    
    t = (1:length(PD_dec))/PD_struct.sampling_freq - PD_struct.basetimes(fo);
    
    for n = 1 %:PD_struct.no_norms
        
        BP_data = load([subj_name, '_wt_BP.mat'], ['BP', PD_struct.norms{n}]);
        
        BP_data = getfield(BP_data, ['BP', PD_struct.norms{n}]);
        
        for ch = 1:2, BP_data(:, :, ch) = zscore(BP_data(:, :, ch)); end
        
        Spec_data = load([subj_name, '_wt.mat'], ['Spec', PD_struct.norms{n}]);
        
        Spec_data = getfield(Spec_data, ['Spec', PD_struct.norms{n}]);
        
        for ty = 2 % 1:PD_struct.no_types
            
            BP_high_cons = getfield(BP_high_consolidated, ['BP', PD_struct.norms{n}, '_high', PD_struct.high_type{ty}]);
            
            for b = 3 % 1:PD_struct.no_bands
                
                %% Get spec & BP data.
                
                interval_length = 10;
                
                BP_hc_length = sum(BP_high_cons(:, b, 3));
                
                BP_hc_blocks = index_to_blocks(BP_high_cons(:, b, 3));
                
                no_blocks = size(BP_hc_blocks, 1);
                
                plot_length = BP_hc_length + (no_blocks - 2)*interval_length;
                
                BP_hc_plot = nan(plot_length, PD_struct.no_bands, 2);
                
                Spec_hc_plot = nan(plot_length, 75, 2);
                
                LFP_hc_plot = nan(plot_length, 2);
                
                t_hc_plot = nan(plot_length, 1);
                
                last_index = 0;
                
                for bl = 1:no_blocks
                    
                    block_length = BP_hc_blocks(bl, 2) - BP_hc_blocks(bl, 1) + 1;
                    
                    block_indices = last_index + (bl > 1)*5 + (1:block_length);
                    
                    t_hc_plot(block_indices) = t(BP_hc_blocks(bl, 1):BP_hc_blocks(bl, 2));
                    
                    BP_hc_plot(block_indices, :, :) = BP_data(BP_hc_blocks(bl, 1):BP_hc_blocks(bl, 2), :, :);
                    
                    Spec_hc_plot(block_indices, :, :) = Spec_data(BP_hc_blocks(bl, 1):BP_hc_blocks(bl, 2), 1:75, :);
                    
                    LFP_hc_plot(block_indices, :) = PD_dec(BP_hc_blocks(bl, 1):BP_hc_blocks(bl, 2), :);
                    
                    last_index = last_index + (bl > 1)*interval_length + block_length;
                    
                end
                
                save([cons_name, '_plot_data.mat'], 't_hc_plot', 'Spec_hc_plot', 'BP_hc_plot', 'LFP_hc_plot')
                
                chunk_plot(cons_name, cons_title, 10, PD_struct, t_hc_plot, Spec_hc_plot, BP_hc_plot, LFP_hc_plot)
                
            end
            
        end
        
    end

end

end

function chunk_plot(save_name, title_name, chunk_size, PD_struct, t, Spec, BP, LFP)

data_length = size(Spec, 1);

no_chunks = ceil(data_length/(chunk_size*PD_struct.sampling_freq));

for c = 1:no_chunks
    
    chunk_name = [save_name, '_chunk', num2str(c, '%03d')];
    
    chunk_start = max(1, (c - 1)*chunk_size*PD_struct.sampling_freq + 1);
    
    chunk_end = min(c*chunk_size*PD_struct.sampling_freq, data_length);
    
    LFP_chunk = LFP(chunk_start:chunk_end, :);
    
    Spec_chunk = Spec(chunk_start:chunk_end, :, :);
    
    BP_chunk = BP(chunk_start:chunk_end, :, :);
    
    t_chunk = t(chunk_start:chunk_end);
    
    t_for_plot = 1:length(t_chunk);
    
    figure
    
    for ch = 1:2
        
        %% Plot Spectrogram.
        
        subplot(4, 1, (ch - 1)*2 + 1)
        
        imagesc(t_for_plot, PD_struct.freqs(1:size(Spec, 2)), Spec_chunk(:, :, ch)')
        
        if ch == 1
            
            title(['Wavelet Spectrogram, ', title_name])
            
        end
        
        set(gca, 'XTickLabel', t_chunk)
        
        ylabel(PD_struct.chan_labels{ch})
        
        axis xy
        
        hold on, plot(t_for_plot', ones(length(t_for_plot), 2)*diag([8 30]), ':w'),
        
        %% Plot LFP & Band Power.
        
        subplot(4, 1, (ch - 1)*2 + 2)
        
        [ax, ~, h2] = plotyy(t_for_plot, LFP_chunk(:, ch), t_for_plot, BP_chunk(:, :, ch));
        
        box off
        
        axis(ax(1), 'tight'), axis(ax(2), 'tight'), xlim(ax(2),xlim(ax(1)))
        
        set(ax(1), 'XTickLabel', t_chunk)
        
        set(ax(2), 'XTickLabel', t_chunk)
        
        set(ax,{'ycolor'},{'black';'black'})
        
        set(get(ax(1),'YLabel'), 'String', 'LFP')
        
        set(get(ax(2),'YLabel'), 'String', 'Band Power')
        
        if ch == 1
            
            legend(h2, PD_struct.band_labels, 'Location', 'SouthWest')
            
        end
        
        axis tight
        
        yl = ylim;
        
        % set(gca, 'YTick', ceil(yl(1)):floor(yl(2)), 'YTickLabel', ceil(yl(1)):2:floor(yl(2)), 'YGrid', 'on')
        
    end
    
    save_as_pdf(gcf, chunk_name)
    
end

end

