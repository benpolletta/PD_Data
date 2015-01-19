function plot_beta_block_datapoints(subject_mat, sd_lim)

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);

load([subject_mat(1:(end - length('_subjects.mat'))), '_', num2str(sd_lim), 'sd_BP_high_dps.mat'])

BP_high_dps(:, :, :, :, :, 2) = BP_high_cum_dps;

for n = 1:no_norms
    
    for t = 1:2
            
        figure
        
        for b = 1:size(BP_high_dps, 2)
            
            for ch = 1:2
                
                subplot(no_bands, 2, 2*b - (2 - ch))
                
                for_plot = reshape(BP_high_dps(:, b, ch, :, n, t), size(BP_high_dps, 1), size(BP_high_dps, 4));
                
                plot([1; 2], for_plot', 's-')
                
                xlim([.5 2.5])
                
            end
            
        end
     
    end
        
end