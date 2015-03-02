function beta_blocks_consolidated_plot(subject_mat, epoch_secs, time_window, percent)

if ~isempty(epoch_secs)
    
    epoch_label = sprintf('_%ddensest', epoch_secs);
    
else
    
    epoch_label = '';
    
end

PD_struct = PD_initialize(subject_mat);

for fo = 1:PD_struct.no_folders
    
    folder = PD_struct.folders{fo};
    
    prefix = PD_struct.prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    cons_name = [subj_name, '_2sd_BP_high_', epoch_label, num2str(time_window/PD_struct.sampling_freq), 's_', num2str(percent), 'pct_consolidated'];
    
    cons_title = [folder, ', High Power, ', num2str(percent), ' Percent of ', num2str(time_window/PD_struct.sampling_freq), ' s'];
    
    load([cons_name, '.mat'], 'BP_high_consolidated')
    
    load([subj_name, '_all_channel_data_dec.mat'])
    
    t = (1:length(PD_dec))/PD_struct.sampling_freq - PD_struct.basetimes(fo);
    
    pd_indices = get_pd_indices(t, fo, subj_name, epoch_secs, BP_suffix, band_index, no_pds);
    
    for n = 1 %:PD_struct.no_norms
        
        BP_data = load([subj_name, '_wt_BP.mat'], ['BP', PD_struct.norms{n}]);
        
        BP_data = getfield(BP_data, ['BP', PD_struct.norms{n}]);
        
        for ch = 1:2, BP_data(:, :, ch) = zscore(BP_data(:, :, ch)); end
        
        Spec_data = load([subj_name, '_wt.mat'], ['Spec', PD_struct.norms{n}]);
        
        Spec_data = getfield(Spec_data, ['Spec', PD_struct.norms{n}]);
        
        for ty = 2 % 1:PD_struct.no_types
            
            BP_high_cons = getfield(BP_high_consolidated, ['BP', PD_struct.norms{n}, '_high', PD_struct.high_type{ty}]);
            
            for b = 3 % 1:PD_struct.no_bands
                
                for pd = 1:PD_struct.no_pds
                    
                    %% Get spec & BP data.
                    
                    interval_length = 10;
                    
                    BP_hc_length = sum(BP_high_cons(:, b, 3));
                    
                    BP_hc_blocks = index_to_blocks(BP_high_cons(:, b, 3) & pd_indices(:, p));
                    
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
                    
                    save([cons_name, PD_struct.pd_label{pd} '_plot_data.mat'], 't_hc_plot', 'Spec_hc_plot', 'BP_hc_plot', 'LFP_hc_plot')
                    
                    chunk_plot(cons_name, cons_title, 10, PD_struct, t_hc_plot, Spec_hc_plot, BP_hc_plot, LFP_hc_plot)
                    
                end
                
            end
            
        end
        
    end

end

end

function chunk_plot(save_name, title_name, chunk_size, PD_struct, t, Spec, BP, LFP)
    
tick_indices = 1:(chunk_size/10):chunk_size;

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
        
        set(ax(1), 'XTick', t_for_plot(tick_indices), 'XTickLabel', t_chunk(tick_indices))
        
        set(ax(2), 'XTick', t_for_plot(tick_indices), 'XTickLabel', t_chunk(tick_indices))
        
        set(ax,{'ycolor'},{'black';'black'})
        
        set(get(ax(1),'YLabel'), 'String', 'LFP')
        
        set(get(ax(2),'YLabel'), 'String', 'Band Power')
        
        if ch == 1
            
            legend(h2, PD_struct.band_labels, 'Location', 'SouthWest')
            
        end
        
        % axis tight
        
        % yl = ylim;
        
        % set(gca, 'YTick', ceil(yl(1)):floor(yl(2)), 'YTickLabel', ceil(yl(1)):2:floor(yl(2)), 'YGrid', 'on')
        
    end
    
    save_as_pdf(gcf, chunk_name)
    
end

end

function pd_indices = get_pd_indices(t, fo, subj_name, epoch_secs, BP_suffix, band_index, no_pds)

if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
    load([subj_name, '_wav_laser_artifacts.mat'], 'laser_periods')
    
    pd_indices = laser_periods;
    
elseif ~isempty(epoch_secs)
    
    if ~isempty(dir([subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs_by_STR.mat']))
        
        load([subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs_by_STR.mat'])
        
        pd_indices = zeros(length(t), no_pds);
        
        for pd = 1:no_pds
            
            bp_max_start = All_bp_max_start(fo, 1, band_index, pd);
            
            bp_max_end = All_bp_max_end(fo, 1, band_index, pd);
            
            pd_indices(bp_max_start:bp_max_end, pd) = 1;
            
        end
        
    elseif ~isempty(dir([subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs.mat']))
        
        load([subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs.mat'])
        
        pd_indices = zeros(length(t), no_pds);
        
        for pd = 1:no_pds
            
            bp_max_start = All_bp_max_start(fo, 1, band_index, pd);
            
            bp_max_end = All_bp_max_end(fo, 1, band_index, pd);
            
            pd_indices(bp_max_start:bp_max_end, pd) = 1;
            
        end
        
    end
    
else
    
    if no_pds == 2
        
        pd_indices(:, 1) = t < 0;
        
        pd_indices(:, 2) = t > 0;
        
    elseif no_pds == 1
        
        pd_indices = ones(length(t), no_pds);
        
    end
    
end

end

