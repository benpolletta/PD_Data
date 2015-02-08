function PD_beta_blocks_rel_infusion_laser_stats(subject_mat, measure, norm_for_power, no_trials_analyzed, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
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

pd_indices = cell(no_folders, 1);

trials = cell(no_folders, no_pds);

no_trials = nan(no_folders, no_pds);

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_wav_laser_artifacts.mat'])
    
    pd_indices{fo} = logical(laser_periods);
    
    for pd = 1:no_pds
        
        trials{fo, pd} = index_to_blocks(laser_periods(:, pd));
        
        no_trials(fo, pd) = size(trials, 1);
        
    end
    
end

if isempty(no_trials_analyzed)

    max_no_trials = all_dimensions(@max, no_trials);

else
    
    max_no_trials = no_trials_analyzed;
    
end
    
no_secs = max_no_trials*5;

% high_type = {'', '_cum'}; no_types = length(high_type);

pct_bp_high = nan(no_secs, no_folders, no_pds, no_bands, no_chans);

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    if strcmp(measure, '_power')
        
        BP_high_cum = get_BP([subj_name, BP_suffix], outlier_lims(fo), norm_for_power, freqs, no_cycles, bands);
        
    elseif isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
        load([subj_name, BP_suffix, '_2sd_BP_high.mat'], 'BP_high_cum')
    
    else
        
        load([subj_name, BP_suffix, '_2sd_BP_high.mat'], 'BP_high')
        
        BP_high_cum = BP_high;
        
    end
    
    length_plotted = min(max_no_trials*triallength(fo)*sampling_freq, size(BP_high_cum, 1));
        
    t = (1:length_plotted)/sampling_freq;
    
    for ch = 1:no_chans
        
        for b = 1:no_bands
            
            figure(b)
            
            beta_blocks_plot = nanconv(BP_high_cum(1:length_plotted, b, ch), ones(sampling_freq, 1)/sampling_freq, 'same');
            
            subplot(no_folders, 2, (fo - 1)*2 + ch)
            
            plot(t/60, beta_blocks_plot)
            
            axis tight
            
            hold on
            
            handle = nan(no_pds, 1);
            
            for pd = 1:no_pds
                
                handle(pd) = plot(t(pd_indices{fo}(1:length_plotted, pd))/60, beta_blocks_plot(pd_indices{fo}(1:length_plotted, pd)), [pd_colors{pd}, '.']);
                
                for tr = 1:min(max_no_trials, no_trials(fo, pd))
                    
                    BP_trial = BP_high_cum(trials{fo, pd}(tr, 1):trials{fo, pd}(tr, 2), b, ch);
                    
                    BP_trial = reshape(BP_trial(1:5*sampling_freq), sampling_freq, 5);
                    
                    pct_bp_high((tr - 1)*5 + (1:5), fo, pd, b, ch) = nansum(BP_trial)'/sampling_freq;
                    
                end
                
            end
            
            axis tight
            
            if fo == 1
                
                title(sprintf('%s, %d - %d Hz', chan_labels{ch}, bands(b, :)))
               
            elseif fo == no_folders
                
                xlabel('Time (Min.)')
                
            end
                    
            if ch == 1
            
                ylabel({folder; 'High Power Density per Sec.'})
                
                if fo == 1
                
                    legend(handle, pd_labels)
                
                end
                    
            end
            
        end
        
    end
    
end

for b = 1:no_bands
           
    save_as_pdf(b, [subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_laser_', num2str(no_trials_analyzed), 'trials_', short_band_labels{b}, measure])
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_laser_', num2str(no_trials_analyzed), 'trials', measure, '.mat'], 'pct_bp_high')

end
