function beta_blocks_rel_infusion_freq_plot(subject_mat, norm, h_norm)
    
close('all')

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

freqs = 1:200;

cycle_lengths = sampling_freq./freqs;

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_pds = length(pd_labels);

no_chans = length(chan_labels);

% pd_labels = {'pre', 'post'};
% 
% pd_colors = {'g', 'r'}; pd_cmap = [0 .5 0; 1 0 0];

fid = nan(no_pds, no_chans);

%% Group average plot (by datapoint). 

[Freq_hist, Spec_high_beta_mean, Spec_high_beta_std] = deal(nan(sum(band_indices{3}), no_pds, no_chans));
    
for ch = 1:no_chans
    
    for pd = 1:no_pds
        
        Freq_data = load([subject_mat(1:(end - length('_subjects.mat'))), '_beta_block_freqs', norm, '_ch', num2str(ch), '_', pd_labels{pd}, '.txt']);
        
        Freq_high_beta = Freq_data(:, 1);
        
        [h, ~] = hist(Freq_high_beta, freqs(band_indices{3}));
        
        if strcmp(h_norm, '_cycles')
            
            Freq_hist(:, pd, ch) = h./cycle_lengths(band_indices{3});
            
        else
        
            Freq_hist(:, pd, ch) = h/sum(h);
        
        end
            
        figure(1)
        
        subplot(1, no_chans, ch)
        
        plot(freqs(band_indices{3}), h/sum(h), pd_colors{pd}, 'LineWidth', 2)
        
        axis tight
        
        xlabel('Freq. (Hz)')
        
        ylabel('Proportion High Beta Datapoints Observed')
        
        if ch == 1
            
            legend(pd_labels)
            
        end
        
        title([chan_labels{ch}, ', Histogram of High Beta Frequencies'])
        
        hold on
        
        Spec_high_beta_mean(:, pd, ch) = nanmean(Freq_data(:, 2:end))';
        
        Spec_high_beta_std(:, pd, ch) = nanstd(Freq_data(:, 2:end))'/sqrt(size(Spec_high_beta_std, 1)/(sampling_freq/8));
        
    end
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), '_beta_block_freqs', norm, h_norm, '_group_stats'], 'Freq_hist', 'Spec_high_beta_mean', 'Spec_high_beta_std')

save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_beta_block_freqs', norm, h_norm, '_group_hist'])

figure

for ch = 1:no_chans
    
    subplot(1, no_chans, ch)

    boundedline(freqs(band_indices{3}), Spec_high_beta_mean(:, :, ch), prep_for_boundedline(Spec_high_beta_std(:, :, ch)), 'cmap', pd_cmap)
      
    axis tight

    xlabel('Freq. (Hz)')
    
    ylabel('Mean \pm S.D. Power')
    
    if ch == 1
        
        legend(pd_labels)
        
    end
    
    title([chan_labels{ch}, ', Mean Power During High Beta'])
    
end

save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_beta_block_freqs', norm, h_norm, '_group_mean'])

%% Group average plot (by recording).

[Freq_hist, Spec_high_beta_mean, Spec_high_beta_std] = deal(nan(no_folders, sum(band_indices{3}), no_pds, no_chans));

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_beta_block_freqs', norm, '.mat'])
    
    if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
        
        load([subj_name, '_wav_laser_artifacts.mat'], 'laser_periods')
        
        pd_indices = laser_periods;
        
    else
        
        pd_indices = nan(length(t), 2);
        
        pd_indices(:, 1) = t < 0;
        
        pd_indices(:, 2) = t > 0;
        
    end
    
    pd_indices = logical(pd_indices);
    
    for ch = 1:no_chans
        
        for pd = 1:no_pds
            
            [h, ~] = hist(Freqs_high_beta(pd_indices(:, pd), ch), freqs(band_indices{3}));
            
            Freq_hist(fo, :, pd, ch) = h'/sum(h);
            
            Spec_high_beta_mean(fo, :, pd, ch) = nanmean(Spec_high_beta(pd_indices(:, pd), :, ch));
            
            Spec_high_beta_std(fo, :, pd, ch) = nanstd(Spec_high_beta(pd_indices(:, pd), :, ch));
            
        end
        
    end
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), '_beta_block_freqs', norm, h_norm, '_individual_stats'], 'Freq_hist', 'Spec_high_beta_mean', 'Spec_high_beta_std')

Freq_hist_mean = permute(mean(Freq_hist), [2 3 4 1]);

Freq_hist_std = permute(std(Freq_hist), [2 3 4 1])/sqrt(no_folders);

figure

for ch = 1:no_chans
    
    subplot(1, no_chans, ch)
    
    boundedline(freqs(band_indices{3})', Freq_hist_mean(:, :, ch),...
        prep_for_boundedline(Freq_hist_std(:, :, ch)), 'cmap', pd_cmap)

    axis tight
    
    xlabel('Freq. (Hz)')
    
    xlim([8 30])
    
    ylabel('Proportion of High Beta')
    
    if ch == 1
        
        legend(pd_labels)
        
    end
    
    title([chan_labels{ch}, ', Histogram of High Beta Frequency'])
   
end

save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_beta_block_freqs', norm, h_norm, '_individual_hist'])

Spec_high_beta_mean_plot = permute(mean(Spec_high_beta_mean), [2 3 4 1]);

Spec_high_beta_std_plot = permute(std(Spec_high_beta_mean), [2 3 4 1])/sqrt(no_folders);

figure

for ch = 1:no_chans
    
    subplot(1, no_chans, ch)
    
    boundedline(freqs(band_indices{3})', Spec_high_beta_mean_plot(:, :, ch),...
        prep_for_boundedline(Spec_high_beta_std_plot(:, :, ch)), 'cmap', pd_cmap)
    
    xlim([8 30])

    axis tight
    
    xlabel('Freq. (Hz)')
    
    ylabel('Mean \pm S.D. Power')
    
    if ch == 1
        
        legend(pd_labels)
        
    end
    
    title([chan_labels{ch}, ', Mean Power During High Beta'])
   
end

save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_beta_block_freqs', norm, h_norm, '_individual_mean'])