function PD_beta_blocks_rel_infusion_pre_post_stats(subject_mat, peak_suffix, epoch_secs, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = peak_suffix;
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end
    
close('all')

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(chan_labels);

no_pds = length(pd_labels);

% pd_colors = {'g', 'r'};

% norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);
% 
% long_norms = {'', ', Increase Over Baseline Power', ', % Total Power', ', Increase in % Total Power Over Baseline'};

% high_type = {'', '_cum'}; no_types = length(high_type);

% if isempty(dir([subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs.mat']))

pct_bp_high = nan(epoch_secs, no_folders, 2, no_bands, 2);

[All_bp_max_start, All_bp_max_end] = deal(nan(no_folders, no_chans, no_bands, no_pds));

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
        load([subj_name, BP_suffix, '_2sd_BP_high.mat'], 'BP_high_cum')
    
    else
        
        load([subj_name, BP_suffix, '_2sd_BP_high.mat'], 'BP_high')
        
        BP_high_cum = BP_high;
        
    end
    
    t = (1:size(BP_high_cum, 1))/sampling_freq - basetimes(fo); % t = t/60;
    
    if no_pds == 2
        
        pd_indices = nan(length(t), 2);
        
        pd_indices(:, 1) = t > (min(t) + epoch_secs/2 - 1) & t < (-epoch_secs/2 + 1);
        
        pd_indices(:, 2) = t > (epoch_secs/2 - 1) & t < (1800 - epoch_secs/2 + 1);
        
    elseif no_pds == 1
        
        pd_indices = ones(length(t), 1);
        
    end
    
    beta_blocks_find = nan(size(BP_high_cum));
    
    for ch = 1:no_chans
        
        for b = 1:no_bands
            
            figure(b)
            
            beta_blocks_find = conv(BP_high_cum(:, b, ch), ones(epoch_secs*sampling_freq, 1)/(epoch_secs*sampling_freq), 'same');
            
            beta_blocks_plot = conv(BP_high_cum(:, b, ch), ones(60*sampling_freq, 1)/(60*sampling_freq), 'same');
            
            subplot(no_folders, 2, (fo - 1)*2 + ch)
            
            plot(t/60, beta_blocks_plot)
            
            axis tight
            
            hold on
            
            handle = nan(no_pds, 1);
            
            for pd = 1:no_pds
                
                % first_possible_index = find(pd_indices(:, pd), 1, 'first') + epoch_secs*sampling_freq/2 - 1;
                % 
                % last_possible_index = find(pd_indices(:, pd), 1, 'last') - epoch_secs*sampling_freq/2 + 1;
                
                [~, bp_max_index] = max(beta_blocks_find(logical(pd_indices(:, pd))));
                
                bp_max_start = max(bp_max_index - floor((epoch_secs/2)*sampling_freq) + find(pd_indices(:, pd), 1, 'first') - 1, 1);
                
                bp_max_end = min(bp_max_start + epoch_secs*sampling_freq, length(t));
                
                handle(pd) = plot(t(bp_max_start:bp_max_end)/60, beta_blocks_plot(bp_max_start:bp_max_end), pd_colors{pd}, 'LineWidth', 2);
                
                for sec = 1:epoch_secs
                   
                    sec_start = max(bp_max_start + (sec - 1)*sampling_freq + 1, 1);
                    
                    sec_end = min(max(bp_max_start + sec*sampling_freq, 1), find(pd_indices(:, pd) == 1, 1, 'last'));
                    
                    pct_bp_high(sec, fo, pd, b, ch) = sum(BP_high_cum(sec_start:sec_end, b, ch))/sampling_freq;
                    
                end
                
                All_bp_max_start(fo, ch, b, pd) = bp_max_start;
                
                All_bp_max_end(fo, ch, b, pd) = bp_max_end;
                
            end
            
            if fo == 1
                
                title(sprintf('%s, %d - %d Hz', chan_labels{ch}, bands(b, :)))
               
            elseif fo == no_folders
                
                xlabel('Time Rel. Infusion (Min.)')
                
            end
                    
            if ch == 1
            
                ylabel({folder; 'High Power Density per Min.'})
                
                if fo == 1
                
                    legend(handle, {['Peak ', num2str(epoch_secs/60), ' Min., Pre-Infusion'], ['Peak ', num2str(epoch_secs/60), ' Min., Post-Infusion']})
                
                end
                    
            end
            
            plot([0; 0], [min(beta_blocks_plot); max(beta_blocks_plot)], 'k', 'LineWidth', 1)
            
        end
        
    end
    
end

for b = 1:no_bands
           
    save_as_pdf(b, [subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_', short_band_labels{b}])
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs.mat'], 'pct_bp_high', 'All_bp_max_start', 'All_bp_max_end')

% else
%     
% load([subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs.mat'])
%     
% end

% if no_pds > 1
%     
%     no_comparisons = nchoosek(no_pds, 2);
%     
%     comparisons = nchoosek(1:no_pds, 2);
%     
% else
%     
%     no_comparisons = 0;
%     
% end