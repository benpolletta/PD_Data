function beta_blocks_rel_infusion_freqs(subject_mat, norm)

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

freqs = 1:200;

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_pds = length(pd_labels);

no_chans = length(chan_labels);

format = make_format(sum(band_indices{3}) + 1, 'f');

fid = nan(no_pds, no_chans);

for pd = 1:no_pds
    
    for ch = 1:no_chans
    
        fid(pd, ch) = fopen([subject_mat(1:(end - length('_subjects.mat'))), '_beta_block_freqs', norm, '_ch', num2str(ch), '_', pd_labels{pd}, '.txt'], 'w');
        
    end
    
end

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_wt.mat'], 'Spec')
    
    load([subj_name, '_2sd_BP_high.mat'], 'BP_high_cum')
    
    t = (1:size(BP_high_cum, 1))/sampling_freq - basetimes(fo); % t = t/60;
    
    if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
        
        load([subj_name, '_wav_laser_artifacts.mat'], 'laser_periods')
        
        pd_indices = laser_periods;
        
    else
        
        pd_indices = nan(length(t), 2);
        
        pd_indices(:, 1) = t < 0;
        
        pd_indices(:, 2) = t > 0;
        
    end
    
    Spec_beta = Spec(:, band_indices{3}, :);
    
    Spec_high_beta = nan(size(Spec_beta));
    
    Freqs_high_beta = nan(size(Spec_beta, 1), no_chans);
        
    for ch = 1:no_chans
        
        Spec_high_beta(logical(BP_high_cum(:, 3, ch)), :, ch) = Spec_beta(logical(BP_high_cum(:, 3, ch)), :, ch);
        
        if strcmp(norm, '_zscore')
            
            [~, Freqs_high_beta(logical(BP_high_cum(:, 3, ch)), ch)] = max(zscore(Spec_high_beta(logical(BP_high_cum(:, 3, ch))), :, ch), [], 2);
            
        else
            
            [~, Freqs_high_beta(logical(BP_high_cum(:, 3, ch)), ch)] = max(Spec_high_beta(logical(BP_high_cum(:, 3, ch)), :, ch), [], 2);
            
        end
        
        Freqs_high_beta(:, ch) = Freqs_high_beta(:, ch) + min(freqs(band_indices{3})) - 1;
        
        for pd = 1:no_pds
            
            fprintf(fid(pd, ch), format, [Freqs_high_beta(logical(pd_indices(:, pd)), ch) Spec_high_beta(logical(pd_indices(:, pd)), :, ch)]');
            
        end
        
    end
    
    save([subj_name, '_beta_block_freqs', norm], 't', 'Spec_high_beta', 'Freqs_high_beta') %, 'Freqs_for_plot')
    
end

fclose('all');