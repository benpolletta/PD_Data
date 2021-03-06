function PD_rel_infusion_plot_spectrogram(subject_mat, peak_suffix, sd_lim, chunk_size, freqs, no_cycles, bands) % , outlier_lims)

% INPUTS:
% subject_mat - .mat file containing info about recordings to be processed.
% sd_lim - standard deviation cutoff above which power in a given frequency
% band is considered "high".
% chunk_size - length of data plotted in each figure (20 sec. is standard).

if isempty(sd_lim), sd_lim = 2; end

if isempty(chunk_size), chunk_size = 20; end

load(subject_mat)

%sampling_freq = 1000;

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200; in_freqs = [];
    
    no_cycles = linspace(3, 21, 200); in_no_cycles = [];
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; in_bands = [];
    
    BP_suffix = ['', peak_suffix];
    
else
    
    in_freqs = freqs; in_no_cycles = no_cycles; in_bands = bands;
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end

no_freqs = length(freqs); no_bands = size(bands, 1);

[band_indices, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands

    % band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    band_labels{b} = [num2str(bands(b, 1)), '-', num2str(bands(b, 2)),'Hz'];

end
    
% BP_chunk = nan(chunk_size*sampling_freq, no_bands); %[total_chunk, beta_chunk] = deal(nan(chunk_size*sampling_freq, 2));

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};  
    
    basetime = basetimes(fo);
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    PD_zs = zscore(PD_dec);
    
    data_length = length(PD_dec);
    
    no_mins = floor(data_length/(sampling_freq*60)); % min(floor(data_length/(sampling_freq*60)), (basetime + 1500)/60);
    
    no_chunks = no_mins*60/chunk_size;
    
    [BP_no_spikes, Spec] = get_BP(subj_name, peak_suffix, [], outlier_lims(fo), '', in_freqs, in_no_cycles, in_bands); % load([subj_name, '_wt.mat'], 'Spec') % 
    
    [BP_norm_no_spikes, ~] = get_BP(subj_name, peak_suffix, [], outlier_lims(fo), '_norm', in_freqs, in_no_cycles, in_bands); % 
    
    load([subj_name, BP_suffix, '_wt_BP.mat'], 'BP', 'BP_norm') % 
    
    load([subj_name, BP_suffix, '_', num2str(sd_lim), 'sd_BP_norm_high.mat'], 'BP_high')
    
    load([subj_name, BP_suffix, '_', num2str(sd_lim), 'sd_BP_high.mat'], 'BP_high_cum')
    
    BP_high(BP_high == 0) = nan; BP_high_cum(BP_high_cum == 0) = nan;
    
    % [BP_zs, BP_norm_zs] = deal(nan(size(BP)));
    
    for ch = 1:2
        
        BP_mean = ones(size(BP(:,:,ch)))*diag(nanmean(BP(:,:,ch)));
        
        BP_std = ones(size(BP(:,:,ch)))*diag(nanstd(BP(:,:,ch)));
        
        BP_no_spikes(:, :, ch) = (BP_no_spikes(:, :, ch) - BP_mean)./BP_std;
    
        BP(:, :, ch) = zscore(BP(:, :, ch));
        
        BP_norm_mean = ones(size(BP_norm(:,:,ch)))*diag(nanmean(BP_norm(:,:,ch)));
        
        BP_norm_std = ones(size(BP_norm(:,:,ch)))*diag(nanstd(BP_norm(:,:,ch)));
        
        BP_norm_no_spikes(:, :, ch) = (BP_norm_no_spikes(:, :, ch) - BP_norm_mean)./BP_norm_std;
        
        BP_norm(:, :, ch) = zscore(BP_norm(:, :, ch));
        
    end
    
    Spike_indicator = peak_loader(subj_name, peak_suffix, size(BP, 1));
        
    for c = 1:no_chunks %(19*3):(29*3) %1:(19*3 - 1) 
        
        chunk_name = [subj_name, BP_suffix, '_chunk', num2str(c, '%03d'), '_wt'];
    
        chunk_start = max(1, (c - 1)*chunk_size*sampling_freq + 1);
        
        chunk_end = min(c*chunk_size*sampling_freq, data_length);
        
        LFP_chunk = PD_zs(chunk_start:chunk_end, :);
        
        WT_chunk = Spec(chunk_start:chunk_end, :, :); % WT_chunk = abs(wavelet_spectrogram(LFP_chunk, sampling_freq, freqs, no_cycles));
        
        BP_chunk = BP(chunk_start:chunk_end, :, :);
        
        BP_norm_chunk = BP_norm(chunk_start:chunk_end, :, :);
        
        Spikes_chunk = Spike_indicator(chunk_start:chunk_end, :);
        
        BP_no_spikes_chunk = BP_no_spikes(chunk_start:chunk_end, :, :);
        
        BP_norm_no_spikes_chunk = BP_norm_no_spikes(chunk_start:chunk_end, :, :);
        
        BP_high_chunk = BP_high(chunk_start:chunk_end, :, :);
        
        BP_high_cum_chunk = BP_high_cum(chunk_start:chunk_end, :, :);
        
        t = ((chunk_start:chunk_end)/sampling_freq - basetime);
        
        minute = ceil(mean(t)/60);
        
        figure
        
        for ch = 1:2
            
            %% Plot Spectrogram.
            
            subplot(6, 1, (ch - 1)*3 + 1)
            
            imagesc(t, freqs(1:100), abs(WT_chunk(:, 1:100, ch))')
            
            if ch == 1
               
                title(['Wavelet Spectrogram, ', folder, ', Minute ', num2str(minute)])
                
            end
            
            ylabel(chan_labels{ch})
            
            axis xy
            
            hold on, plot(t', ones(length(t), 2)*diag([10 30]), ':w'),
            
            % subplot(6, 1, (ch - 1)*3 + 2)
            % 
            % imagesc(t, freqs(1:100), zscore(WT_chunk(:, 1:100, ch))')
            % 
            % ylabel('Freq. Normalized')
            % 
            % axis xy
            % 
            % hold on, plot(t', ones(length(t), 2)*diag([10 30]), '--w'),
            % 
            % % hold on, plot(t, 100*(LFP_chunk(:, ch) - min(LFP_chunk(:, ch)))/range(LFP_chunk(:, ch)), '--r')
            
            %% Plot LFP & Band Power.
            
            subplot(6, 1, (ch - 1)*3 + 2)
            
            [ax, ~, h2] = plotyy(t, LFP_chunk(:, ch), t, BP_no_spikes_chunk(:, :, ch));
            
            box off
            
            hold(ax(1), 'on'), hold(ax(2), 'on')
            
            plot(ax(1), t(logical(Spikes_chunk(:, ch))), LFP_chunk(logical(Spikes_chunk(:, ch)), ch), 'v')
            
            plot(ax(2), t, BP_chunk(:, :, ch), ':')
            
            plot(ax(2), t, BP_high_cum_chunk(:, :, ch).*BP_chunk(:, :, ch), 'LineWidth', 2)
            
            axis(ax(1), 'tight'), axis(ax(2), 'tight'), xlim(ax(2),xlim(ax(1)))
            
            set(ax,{'ycolor'},{'black';'black'})
            
            set(get(ax(1),'YLabel'), 'String', 'LFP')
            
            set(get(ax(2),'YLabel'), 'String', 'Band Power')
            
            if ch == 1
            
                legend(h2, band_labels, 'Location', 'SouthWest')

            end
                
            axis tight
            
            yl = ylim;
            
            set(gca, 'YTick', ceil(yl(1)):floor(yl(2)), 'YTickLabel', ceil(yl(1)):2:floor(yl(2)), 'YGrid', 'on')
            
            % %% Plot LFP & Band Power.
            % 
            % subplot(6, 1, (ch - 1)*3 + 3)
            % 
            % plot(t, [LFP_chunk(:, ch) BP_zs(chunk_start:chunk_end, :, ch)])%BP_chunk]))
            % 
            % axis tight
            % 
            % yl = ylim;
            % 
            % set(gca, 'YTick', ceil(yl(1)):floor(yl(2)), 'YTickLabel', ceil(yl(1)):2:floor(yl(2)), 'YGrid', 'on')
            % 
            % ylabel(['LFP & Band Power'])
            % 
            % if ch == 1
            % 
            %     legend({'LFP', band_labels{:}})
            % 
            % end
            
            %% Plot LFP & Band Power (Normalized by Total Power).
            
            subplot(6, 1, (ch - 1)*3 + 3)
            
            [ax, ~, ~] = plotyy(t, nan(size(LFP_chunk(:, ch))), t, BP_norm_no_spikes_chunk(:, :, ch));
            
            box off
            
            hold(ax(1), 'on'), hold(ax(2), 'on')
            
            plot(ax(2), t, BP_norm_chunk(:, :, ch), ':')
            
            plot(ax(2), t, BP_high_chunk(:, :, ch).*BP_norm_chunk(:, :, ch), 'LineWidth', 2)
            
            % plot(ax(2), t, BP_chunk(:, :, ch), ':')
            % 
            % plot(ax(2), t, BP_high_chunk(:, :, ch).*BP_chunk(:, :, ch), 'LineWidth', 2)
            
            axis(ax(1), 'tight'), axis(ax(2), 'tight'), xlim(ax(2),xlim(ax(1)))
            
            set(ax,{'ycolor'},{'black';'black'})
            
            set(get(ax(1),'YLabel'), 'String', 'LFP')
            
            set(get(ax(2),'YLabel'), 'String', '% Band Power')
                
            axis tight
            
            yl = ylim;
            
            set(gca, 'YTick', ceil(yl(1)):floor(yl(2)), 'YTickLabel', ceil(yl(1)):2:floor(yl(2)), 'YGrid', 'on')
            
        end
        
        % %% Plot LFP.
        % 
        % subplot(5, 1, 5)
        % 
        % [ax, h1, h2] = plotyy(t, LFP_chunk, t, beta_chunk./total_chunk);
        % 
        % %title(['Raw LFP & Pct. \beta Power, ', folder, ', Minute ', num2str(minute)])
        % 
        % box off
        % 
        % axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1))) %axis tight
        % 
        % set(ax,{'ycolor'},{'black';'black'})
        % 
        % set(get(ax(1),'YLabel'),'String','LFP')
        % 
        % set(get(ax(2),'YLabel'),'String','% \beta Power')
        % 
        % legend(h1, chan_labels, 'Location', 'NorthWest')
        % 
        % legend(h2, chan_labels, 'Location', 'SouthWest')
        
        save_as_pdf(gcf, chunk_name)
        
    end
    
    close('all')
    
    % end
    
end

end