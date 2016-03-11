function collect_control_opto_freqs(time_window, percent, norm, band_index, freqs, no_cycles, bands)

sampling_freq = 500;

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

subject_matnames = {'st_m1_ali_control_opto', 'st_m1_ali2_control_opto', 'st_m1_ali3_control_opto'};

no_mats = length(subject_matnames);

channels = [1 2; 2 1; 2 1]; no_chans = size(channels, 2);

load([subject_matnames{1}, '_subjects.mat'], 'pd_labels')

no_pds = length(pd_labels); % no_chans = length(chan_labels);

format = make_format(sum(band_indices{band_index}) + 2, 'f');

fid = nan(no_pds, no_chans);

for pd = 1:no_pds
    
    for ch = 1:no_chans
    
        fid(pd, ch) = fopen(['CONTROL_OPTO', BP_suffix, '_2sd_', short_band_labels{band_index}, '_high_',...
        num2str(time_window/sampling_freq), 's_', num2str(percent), 'pct_consolidated_freqs',...
        norm, '_ch', num2str(ch), '_', pd_labels{pd}, '.txt'], 'w');
        
    end
    
end

for s = 1:no_mats
    
    for pd = 1:no_pds
        
        for ch = 1:no_chans
            
            subj_pd_ch_freqs = load([subject_matnames{s}, BP_suffix, '_2sd_', short_band_labels{band_index}, '_high_',...
                num2str(time_window/sampling_freq), 's_', num2str(percent), 'pct_consolidated_freqs',...
                norm, '_ch', num2str(channels(s, ch)), '_', pd_labels{pd}, '.txt']);
            
            fprintf(fid(pd, ch), format, subj_pd_ch_freqs');
    
        end
        
    end
    
end

fclose('all');