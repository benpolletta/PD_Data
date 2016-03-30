function PD_beta_blocks_rel_infusion_pre_post_power(subject_mat, peak_suffix, epoch_secs, pd_handle, norm, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200; in_freqs = [];
    
    no_cycles = linspace(3, 21, length(freqs)); in_no_cycles = [];
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; in_bands = [];
    
    BP_suffix = ['', peak_suffix];
    
else
    
    in_freqs = freqs; in_no_cycles = no_cycles; in_bands = bands;
    
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

BP_sec = nan(epoch_secs, no_folders, 2, no_bands, 2);

load([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, '.mat'])

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    outlier_lim = outlier_lims(fo);
    
    subj_name = [folder,'/',prefix];
    
    if strcmp(norm, '_peaks')
        
        data = load([subj_name, BP_suffix, '_wt_BP.mat']);
        
        Spike_indicator = peak_loader(subj_name, peak_suffix, size(data.BP, 1));
        
        BP = repmat(permute(Spike_indicator, [1 2 3]), [1 no_bands 1]);
        
    else
        
        BP = get_BP(subj_name, peak_suffix, outlier_lim, norm, in_freqs, in_no_cycles, in_bands);
    
        % load([subj_name, BP_suffix, '_wt_BP.mat'])
        % 
        % eval(['BP = BP', norm, ';'])
        
    end
    
    t = (1:size(BP, 1))/sampling_freq - basetimes(fo); % t = t/60;
    
    clear pd_indices
    
    if no_pds == 2
        
        pd_indices(:, 1) = t < 0; pd_indices(:, 2) = t > 0;
        
    elseif no_pds == 1
        
        pd_indices = ones(length(t), 1);
        
    end
    
    pd_indices = logical(pd_indices);
    
    % if ~strcmp(norm, '_peaks')
    % 
    %     if ~isempty(outlier_lims) && ~isempty(dir([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat']))
    % 
    %         load([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat'])
    % 
    %         [~, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, no_cycles, bands);
    % 
    %         outlier_nans = repmat(outlier_nans, [1 1 2]);
    % 
    %         BP(logical(outlier_nans)) = nan;
    % 
    %     end
    % 
    %     Spike_indicator = peak_loader(subj_name, peak_suffix, size(BP, 1));
    % 
    %     [~, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
    % 
    %     BP(logical(spike_nans)) = nan;
    % 
    % end
    
    % beta_blocks_find = nan(size(BP_high_cum));
    
    for b = 1:no_bands
        
        figure(b)
        
        for ch = 1:no_chans
            
            handle = nan(no_pds, 1);
            
            for pd = 1:no_pds
                
                bp_max_start = All_bp_max_start(fo, ch, b, pd);
                
                bp_max_end = All_bp_max_end(fo, ch, b, pd);
                
                BP_plot = nanconv(BP(:, b, ch), ones(60*sampling_freq, 1)/(60*sampling_freq), 'nanout');
                
                subplot(no_folders, 2, (fo - 1)*2 + ch)
                
                plot(t/60, BP_plot)
                
                axis tight
                
                hold on
                
                handle(pd) = plot(t(bp_max_start:bp_max_end)/60, BP_plot(bp_max_start:bp_max_end), pd_colors{pd}, 'LineWidth', 2);
                
                BP_selected = BP(bp_max_start:bp_max_end, b, ch);
                
                BP_selected = nans_to_end(BP_selected);
                
                for sec = 1:epoch_secs
                    
                    sec_start = (sec - 1)*sampling_freq + 1;
                    
                    sec_end = min(sec*sampling_freq, size(BP_selected, 1));
                    
                    % sec_start = max(bp_max_start + (sec - 1)*sampling_freq + 1, 1);
                    % 
                    % sec_end = min(max(bp_max_start + sec*sampling_freq, 1), find(pd_indices(:, pd) == 1, 1, 'last'));
                    
                    BP_sec(sec, fo, pd, b, ch) = nanmean(BP_selected(sec_start:sec_end));
                    
                end
                
            end
            
            if fo == 1
                
                title(sprintf('%s, %d - %d Hz', chan_labels{ch}, bands(b, :)))
                
            elseif fo == no_folders
                
                xlabel('Time Rel. Infusion (Min.)')
                
            end
            
            if ch == 1
                
                ylabel({folder; 'High Power Density per Min.'})
                
                if fo == 1
                    
                    for pd = 1:no_pds
                        
                        ledge{pd} = ['Peak ', num2str(epoch_secs/60), ' Min., ', pd_labels{pd}];
                        
                    end
                    
                    legend(handle, ledge)
                    
                end
                
            end
            
            plot([0; 0], [min(BP_plot); max(BP_plot)], 'k', 'LineWidth', 1)
            
        end
        
    end
    
end

for b = 1:no_bands
           
    save_as_pdf(b, [subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_', short_band_labels{b}, pd_handle, norm, '_power'])
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, norm, '_power.mat'], 'BP_sec')

end
