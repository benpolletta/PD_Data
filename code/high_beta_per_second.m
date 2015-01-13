function high_beta_per_second(subject_mat, sd_lim)

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

band_labels = cell(no_bands, 1);

for b = 1:no_bands
   
    band_labels{b} = sprintf('%d - %d Hz', bands(b, 1), bands(b, 2));
     
end

norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);

high_type = {'', '_cum'}; no_types = length(high_type);

format = make_format(no_bands + 1, 'f');
    
for n = 1:no_norms
    
    for ty = 1:no_types
        
        for ch = 1:2
            
            fid(ch) = fopen([subject_mat(1:(end - length('_subjects.mat'))), '_ch', num2str(ch), '_BP', norms{n}, '_high', high_type{ty}, '_per_sec.txt'], 'w');
            
        end
        
        for fo = 1:no_folders
            
            folder = folders{fo};
            
            prefix = prefixes{fo};
            
            subj_name = [folder,'/',prefix];
            
            BP_high_data = load([subj_name, '_', num2str(sd_lim), 'sd_BP', norms{n}, '_high.mat'], ['BP_high', high_type{ty}]);
            
            BP_high_data = getfield(BP_high_data, ['BP_high', high_type{ty}]);
            
            no_secs = ceil(length(BP_high_data)/sampling_freq);
            
            BP_high_per_sec = nan(no_secs, no_bands, 2);
            
            BP_high_data((end + 1):(no_secs*sampling_freq), :, :) = 0;
            
            t = (1:size(BP_high_data, 1))/sampling_freq - basetimes(fo);
            
            for ch = 1:2
                
                for b = 1:no_bands
                    
                    BP_high_band_chan = BP_high_data(:, b, ch);
                    
                    BP_high_per_sec(:, b, ch) = nansum(reshape(BP_high_band_chan, sampling_freq, no_secs))'/sampling_freq;
                    
                end
                
                t_sec = nanmean(reshape(t, sampling_freq, no_secs))';
                
                fprintf(fid(ch), format, [t_sec BP_high_per_sec(:, :, ch)]');
                
            end
            
        end
        
        for ch = 1:2, fclose(fid(ch)); end
        
    end
    
end