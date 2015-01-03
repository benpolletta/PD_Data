function PD_beta_blocks_rel_infusion_plot(subject_mat, sd_lim, epoch_secs)
    
load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];

band_labels = cell(length(bands), 1);

for b = 1:length(bands)
   
    band_labels{b} = sprintf('%d - %d Hz', bands(b, 1), bands(b, 2));
     
end

% norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);

% high_type = {'', '_cum'}; no_types = length(high_type);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    % base_index = basetime*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    % BP_data = load([subj_name,'_wt_BP.mat']);
    
    % BP_data = getfield(BP_data, ['BP', norms{n}]);
    %
    % BP_data = zscore(BP_data);
    %
    % t = (1:size(BP_data, 1))/sampling_freq;
    
    % load([subj_name, '_', num2str(sd_lim), 'sd_BP_norm_high.mat'], 'BP_high')
    
    load([subj_name, '_', num2str(sd_lim), 'sd_BP_high.mat'], 'BP_high', 'BP_high_cum')
    
    t = (1:size(BP_high, 1))/sampling_freq - basetimes(fo);
    
    beta_blocks = nan([size(BP_high), 2]);
    
    beta_blocks(:, :, :, 1) = BP_high;
    
    beta_blocks(:, :, :, 2) = BP_high_cum;
    
    beta_blocks_plot = beta_blocks;
    
    no_bands = size(beta_blocks, 2);
    
    for i = 1:2
        
        figure(i)
        
        for ch = 1:2
            
            for b = 1:no_bands
                
                beta_blocks_plot(:, b, ch, i) = zscore(conv(beta_blocks(:, b, ch, i), ones(epoch_secs*sampling_freq, 1), 'same'));
                
            end
            
            subplot(no_folders, 2, 2*fo - (2 - ch))
            
            plot(t, beta_blocks_plot(:, :, ch, i))
            
            axis tight
            
            if ch == 1 
                
                ylabel(folder)
                
                if fo == 1
            
                    legend(band_labels)
                    
                end
            
            end
            
            title(chan_labels{ch})
            
            hold on
            
            plot([basetime; basetime], [all_dimensions(@min, beta_blocks_plot(:, :, ch, i)); all_dimensions(@max, beta_blocks_plot(:, :, ch, i))], 'r')
            
        end
        
    end
    
    % figure
    %
    % for ch = 1:2
    %
    % subplot(2, 1, ch)
    %
    % plot(t, BP_data(:, :, ch))
    
    % t_mat = repmat(t, 1, size(BP_data, 2));
    %
    % hold on
    %
    % plot(t_mat(logical(BP_high(:, :, ch))), BP_data(logical(BP_high(:, :, ch)), :, ch), 'LineWidth', 2)
    %
    % plot(t(BP_high_cum(:, :, ch)), BP_data(BP_high_cum(:, :, ch), :, ch), ':', 'LineWidth', 2)
    %
    % plot(t(BP_high(:, :, ch) & BP_high_cum(:, :, ch)), BP_data(BP_high(:, :, ch) & BP_high_cum(:, :, ch), :, ch), 'LineWidth', 3)
    
end

end