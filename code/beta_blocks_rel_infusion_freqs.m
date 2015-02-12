function beta_blocks_rel_infusion_freqs(subject_mat, norm, band_index, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

load(subject_mat)

subj_mat_name = subject_mat(1:(end - length('_subjects.mat')));

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, BP_suffix, '_wt.mat'], 'sampling_freq')

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_pds = length(pd_labels);

no_chans = length(chan_labels);

format = make_format(sum(band_indices{band_index}) + 1, 'f');

fid = nan(no_pds, no_chans);

for pd = 1:no_pds
    
    for ch = 1:no_chans
    
        fid(pd, ch) = fopen([subj_mat_name, BP_suffix, '_beta_block_', short_band_labels{band_index}, '_freqs', norm, '_ch', num2str(ch), '_', pd_labels{pd}, '.txt'], 'w');
        
    end
    
end

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    if strcmp(norm, '_pct')
        
        load([subj_name, BP_suffix, '_wt_pct.mat'], 'Spec_pct')
        
        Spec = Spec_pct;
        
    else
    
        load([subj_name, BP_suffix, '_wt.mat'], 'Spec')
        
        Spec = abs(Spec);
    
    end
        
    load([subj_name, BP_suffix, '_2sd_BP_high.mat'], 'BP_high_cum')
    
    t = (1:size(BP_high_cum, 1))/sampling_freq - basetimes(fo); % t = t/60;
    
    if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
        
        load([subj_name, '_wav_laser_artifacts.mat'], 'laser_periods')
        
        pd_indices = laser_periods;
        
    elseif ~isempty(dir([subj_mat_name, BP_suffix, '_pct_BP_high_2.5_min_secs_by_STR.mat']))
        
        % load([subj_mat_name, BP_suffix, '_pct_BP_high_2.5_min_secs_by_STR.mat'])
        % 
        % pd_indices = zeros(length(t), no_pds);
        % 
        % for pd = 1:no_pds
        % 
        % bp_max_start = All_bp_max_start(fo, 1, band_index, pd);
        % 
        % bp_max_end = All_bp_max_end(fo, 1, band_index, pd);
        % 
        % bp_max_index = zeros(size(t));
        % 
        % pd_indices(bp_max_start:bp_max_end, pd) = 1;
        % 
        % end
        
        pd_indices(:, 1) = t < 0;
        
        pd_indices(:, 2) = t > 0;
        
    elseif ~isempty(dir([subj_mat_name, BP_suffix, '_pct_BP_high_2.5_min_secs.mat']))
        
        % load([subj_mat_name, BP_suffix, '_pct_BP_high_2.5_min_secs.mat'])
        % 
        % pd_indices = zeros(length(t), no_pds);
        % 
        % for pd = 1:no_pds
        % 
        % bp_max_start = All_bp_max_start(fo, 1, band_index, pd);
        % 
        % bp_max_end = All_bp_max_end(fo, 1, band_index, pd);
        % 
        % bp_max_index = zeros(size(t));
        % 
        % pd_indices(bp_max_start:bp_max_end, pd) = 1;
        % 
        % end
        
        pd_indices = ones(length(t), no_pds);
        
    end
    
    Spec_beta = Spec(:, band_indices{band_index}, :);
    
    Spec_high_beta = nan(size(Spec_beta));
    
    Freqs_high_beta = nan(size(Spec_beta, 1), no_chans);
        
    for ch = 1:no_chans
        
        Spec_high_beta(logical(BP_high_cum(:, band_index, ch)), :, ch) = Spec_beta(logical(BP_high_cum(:, band_index, ch)), :, ch);
        
        if strcmp(norm, '_zscore')
            
            Spec_hb_zs = nanzscore(Spec_high_beta(:, :, ch));
            
            [~, Freqs_high_beta(:, ch)] = nanmax(Spec_hb_zs, [], 2);
            
        else
            
            [~, Freqs_high_beta(:, ch)] = nanmax(Spec_high_beta(:, :, ch), [], 2);
            
        end
        
        Freqs_high_beta(:, ch) = Freqs_high_beta(:, ch) + min(freqs(band_indices{band_index})) - 1;
        
        for pd = 1:no_pds
            
            fprintf(fid(pd, ch), format, [Freqs_high_beta(logical(pd_indices(:, pd)), ch) Spec_high_beta(logical(pd_indices(:, pd)), :, ch)]');
            
        end
        
    end
    
    save([subj_name, BP_suffix, '_beta_block_', short_band_labels{band_index}, '_freqs', norm], 't', 'Spec_high_beta', 'Freqs_high_beta') %, 'Freqs_for_plot')
    
end

fclose('all');