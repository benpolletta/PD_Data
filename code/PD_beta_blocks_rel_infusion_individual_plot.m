function PD_beta_blocks_rel_infusion_individual_plot(folder, prefix, basetime, sd_lim, epoch_secs)

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];

band_labels = cell(length(bands), 1);

for b = 1:length(bands)
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, 1), bands(b, 2));
    
end

high_labels = {'BP >= 2 S.D.', 'BP >= 2 S.D., No Lower High BP'};

subj_name = [folder,'/',prefix];

load([subj_name, '_wt_BP.mat'], 'sampling_freq')

load([subj_name, '_', num2str(sd_lim), 'sd_BP_high.mat'], 'BP_high', 'BP_high_cum')

t = (1:size(BP_high, 1))/sampling_freq - basetime;

beta_blocks = nan([size(BP_high), 2]);

beta_blocks(:, :, :, 1) = BP_high;

beta_blocks(:, :, :, 2) = BP_high_cum;

beta_blocks_plot = beta_blocks;

no_bands = size(beta_blocks, 2);

figure

for i = 1:2
    
    for ch = 1:2
        
        for b = 1:no_bands
            
            beta_blocks_plot(:, b, ch, i) = zscore(conv(beta_blocks(:, b, ch, i), ones(epoch_secs*sampling_freq, 1), 'same'));
            
        end
        
        subplot(2, 2, 2*i - (2 - ch))
        
        plot(t, beta_blocks_plot(:, :, ch, i))
        
        axis tight
        
        if ch == 1
            
            ylabel(folder)
            
            if i == 1
                
                legend(band_labels)
                
            end
            
        end
        
        title(['Channel ', num2str(ch), ', ', high_labels{i}])
        
        hold on
        
        plot([basetime; basetime], [all_dimensions(@min, beta_blocks_plot(:, :, ch, i)); all_dimensions(@max, beta_blocks_plot(:, :, ch, i))], 'r')
        
    end
    
end

end