function high_beta_per_sec_histogram(subject_mat)

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

pd_colors = {'g', 'r'};
    
for n = 1:no_norms
    
    for ty = 1:no_types
        
        figure
        
        for ch = 1:2
            
            BP_high_per_sec_data = load([subject_mat(1:(end - length('_subjects.mat'))), '_ch', num2str(ch), '_BP', norms{n}, '_high', high_type{ty}, '_per_sec.txt']);
            
            t_sec = BP_high_per_sec_data(:, 1);
            
            pd_index = nan(length(t_sec), 2);
                
            pd_index(:, 1) = t_sec < 0; pd_index(:, 2) = t_sec > 0;
            
            BP_high_per_sec = BP_high_per_sec_data(:, 2:end);
            
            for b = 1:no_bands
               
                BP_high_per_sec_band = BP_high_per_sec(:, b);
                
                range = [min(BP_high_per_sec_band) max(BP_high_per_sec_band)];
                
                range = max(eps, range);
                
                [bin_edges, bin_centers] = make_bins(range(1), range(2), 100, 'log');
                
                subplot(no_bands, 2, (b - 1)*2 + ch)
                
                for pd = 1:2
                   
                    [h, ~] = histc(BP_high_per_sec_band(logical(pd_index(:, pd))), bin_edges);
                    
                    % h = max(eps, h);
                    
                    loglog(bin_centers, h(1:(end - 1))/sum(pd_index(:, pd)), pd_colors{pd})
                    
                    hold on
                    
                end
                
                % axis tight
                
            end
        
        end
        
    end
    
end