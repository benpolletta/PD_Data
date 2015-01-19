function PD_beta_blocks_rel_infusion_pre_post_plot_OLD(subject_mat, sd_lim, epoch_secs)
    
load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

band_labels = cell(no_bands, 1);

for b = 1:no_bands
   
    band_labels{b} = sprintf('%d - %d Hz', bands(b, 1), bands(b, 2));
     
end

norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);

long_norms = {'', ', Increase Over Baseline Power', ', % Total Power', ', Increase in % Total Power Over Baseline'};

high_type = {'', '_cum'}; no_types = length(high_type);

[bp_max, bp_min] = deal(nan(no_folders, 2, size(bands, 1), 2, 2));
    
for n = 1:no_norms
    
    for ty = 1:no_types
        
        for fo = 1:no_folders
            
            folder = folders{fo};
            
            prefix = prefixes{fo};
            
            subj_name = [folder,'/',prefix];
            
            load([subj_name, '_', num2str(sd_lim), 'sd_BP', norms{n}, '_high.mat'], 'BP_high', 'BP_high_cum')
            
            t = (1:size(BP_high, 1))/sampling_freq - basetimes(fo);
            
            pd_indices = nan(length(t), 2);
            
            pd_indices(:, 1) = t < 0; pd_indices(:, 2) = t > 0 & t < 1800;
            
            beta_blocks = nan([size(BP_high), 2]);
            
            beta_blocks(:, :, :, 1) = BP_high;
            
            beta_blocks(:, :, :, 2) = BP_high_cum;
            
            beta_blocks_plot = beta_blocks;
            
            % no_bands = size(beta_blocks, 2);
            
            for ch = 1:2
                
                for b = 1:no_bands
                    
                    beta_blocks_plot(:, b, ch, ty) = conv(beta_blocks(:, b, ch, ty), ones(epoch_secs*sampling_freq, 1)/(epoch_secs*sampling_freq), 'same');
                    
                    for pd = 1:2
                        
                        bp_max(fo, pd, b, ch, ty) = max(beta_blocks_plot(logical(pd_indices(:, pd)), b, ch, ty));
                        
                        bp_min(fo, pd, b, ch, ty) = min(beta_blocks_plot(logical(pd_indices(:, pd)), b, ch, ty));
                        
                    end
                    
                end
                
            end
            
        end
        
        figure((n - 1)*2 + ty)
        
        for b = 1:no_bands
            
            for ch = 1:2
                
                subplot(no_bands, 4, (b - 1)*4 + (ch - 1)*2 + 1)
                
                plot([1; 2], bp_max(:, :, b, ch, ty)', 's-')
                
                xlim([.5 2.5])
                
                subplot(no_bands, 4, (b - 1)*4 + (ch - 1)*2 + 2)
                
                plot([1; 2], bp_min(:, :, b, ch, ty)', 's-')
                
                xlim([.5 2.5])
                
            end
            
        end
    
        save_as_pdf((n - 1)*2 + ty, [subject_mat(1:(end - length('_subjects.mat'))), '_wt_BP_high', norms{n}, high_type{ty}, '_pre_post_', num2str(epoch_secs/60), '_mins'])
    
    end
    
end

end