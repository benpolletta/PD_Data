function PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual(subject_mat, peak_suffix, epoch_secs, pd_handle, norm, band_index_for_time, band_index_for_display, freqs, no_cycles, bands)

% Leave epoch_secs empty when using for optogenetics data, and enter
% '_ntrials' for the argument pd_handle.

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = ['', peak_suffix];
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end
    
close('all')

load(subject_mat)

no_folders = length(folders);

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
    
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

display_indices = band_indices{band_index_for_display};

no_chans = length(chan_labels);

no_pds = length(pd_labels);

subj_mat_name = subject_mat(1:(end - length('_subjects.mat')));

if ~isempty(epoch_secs)
    
    spectrum_name = [subj_mat_name, BP_suffix, '_pct_', short_band_labels{band_index_for_time}, '_high_',...
        num2str(epoch_secs/60), '_min_secs', pd_handle, norm, '_spectrum'];
    
else
    
    spectrum_name = [subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_',...
        short_band_labels{band_index_for_time}, pd_handle, norm, '_spectrum'];
    
    epoch_secs = str2double(pd_handle(2:(end - length('trials'))));
    
end

load([spectrum_name, '.mat'])

no_freqs = sum(display_indices); % size(WT_sec, 1);
        
% if strcmp(norm, '_pct')
% 
%     freq_multiplier = ones(no_freqs, size(WT_sec, 2));
% 
% else
% 
%     freq_multiplier = repmat(freqs', 1, size(WT_sec, 2));
% 
% end

%% Individual bar plots.

figure

All_mean = nan(no_freqs, no_folders, no_pds, no_chans);

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder, '/', prefix];
    
    for ch = 1:no_chans
        
        [WT_mean, WT_se] = deal(nan(no_freqs, no_pds));
        
        for pd = 1:2 % no_pds
            
            WT_mean(:, pd) = nanmean(WT_sec(display_indices, :, fo, pd, ch), 2); % .*freq_multiplier, 2);
            
            WT_se(:, pd) = nanstd(WT_sec(display_indices, :, fo, pd, ch), [], 2)/sqrt(epoch_secs); % .*freq_multiplier, [], 2)/sqrt(epoch_secs);
            
        end
        
        subplot(no_chans, no_folders, (ch - 1)*no_folders + fo)
        % subplot(no_folders, no_chans, (fo - 1)*no_chans + ch) % subplot(no_bands, 2, (b - 1)*2 + ch)
        
        if strcmp(norm, '_pct')
            
            boundedline(freqs(display_indices), WT_mean, prep_for_boundedline(norminv(1 - .05, 0, 1)*WT_se))
            
        else
            
            freq_multiplier = repmat(freqs(display_indices)', 1, no_pds);
            
            boundedline(freqs(display_indices), WT_mean.*freq_multiplier, prep_for_boundedline(WT_se.*freq_multiplier))
            
        end
        
        axis tight
        
        % xlim([8 30])
        
        All_mean(:, fo, :, ch) = permute(WT_mean, [1 3 2]);
        
        if fo == 1
            
            title({[num2str(epoch_secs/60), ' Minutes of Densest High Power'];[folder, ', ', band_labels{b}]})
            % title({chan_labels{ch};[num2str(epoch_secs/60), ' Minutes of Densest High Power'];[folder, ', ', band_labels{b}]})
            
            ylabel(chan_labels{ch})
            
            legend(pd_labels(1:2))
            
        else
            
            title(folder)
            % ylabel(folder)
            
        end
        
    end
    
end

save_as_pdf(gcf, [spectrum_name, '_', short_band_labels{band_index_for_display}, '_individual'])

%% Stats & figures treating each individual as an observation.

figure

for ch = 1:no_chans
    
    for pd = 1:no_pds
        
        All_mean_mean(:, pd) = nanmean(All_mean(:, :, pd, ch), 2);
        
        All_mean_se(:, pd) = nanstd(All_mean(:, :, pd, ch), [], 2)/sqrt(no_folders);
        
    end
    
    subplot(1, no_chans, ch)
    
    if strcmp(norm, '_pct')
    
        boundedline(freqs(display_indices), All_mean_mean, prep_for_boundedline(norminv(1 - .05, 0, 1)*All_mean_se))
        
    else
            
        freq_multiplier = repmat(freqs(display_indices)', 1, no_pds);
    
        boundedline(freqs(display_indices), All_mean_mean.*freq_multiplier, prep_for_boundedline(All_mean_se.*freq_multiplier))
    
    end
        
    axis tight
    
    title([chan_labels{ch}, ', ', num2str(epoch_secs/60), ' Minutes of Densest High Power, ', band_labels{band_index_for_time}])
    
    legend(pd_labels(1:2))
    
end

save_as_pdf(gcf, [spectrum_name, '_', short_band_labels{band_index_for_display}, '_individual_avg'])

end
