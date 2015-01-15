function PD_beta_blocks_rel_infusion_plot(subject_mat, sd_lim, epoch_secs)
    
load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];

band_labels = cell(length(bands), 1);

for b = 1:length(bands)
   
    band_labels{b} = sprintf('%d - %d Hz', bands(b, 1), bands(b, 2));
     
end

no_norms = length(norms);

high_type = {'', '_cum'}; no_types = length(high_type);
    
for n = 1:no_norms
    
    for ty = 1:no_types
        
        figure((n - 1)*2 + ty)
        
        for fo = 1:length(folders)
            
            folder = folders{fo};
            
            prefix = prefixes{fo};
            
            % basetime = basetimes(fo);
            
            % base_index = basetime*sampling_freq;
            
            subj_name = [folder,'/',prefix];
            
            % BP_data = load([subj_name,'_wt_BP.mat']);
            
            % BP_data = getfield(BP_data, ['BP', norms{n}]);
            %
            % BP_data = zscore(BP_data);
            %
            % t = (1:size(BP_data, 1))/sampling_freq;
            
            % load([subj_name, '_', num2str(sd_lim), 'sd_BP_norm_high.mat'], 'BP_high')
            
            load([subj_name, '_', num2str(sd_lim), 'sd_BP', norms{n}, '_high.mat'], 'BP_high', 'BP_high_cum')
            
            t = (1:size(BP_high, 1))/sampling_freq - basetimes(fo);
            
            beta_blocks = nan([size(BP_high), 2]);
            
            beta_blocks(:, :, :, 1) = BP_high;
            
            beta_blocks(:, :, :, 2) = BP_high_cum;
            
            beta_blocks_plot = beta_blocks;
            
            no_bands = size(beta_blocks, 2);
            
            for ch = 1:2
                
                for b = 1:no_bands
                    
                    beta_blocks_plot(:, b, ch, ty) = zscore(conv(beta_blocks(:, b, ch, ty), ones(epoch_secs*sampling_freq, 1)/(epoch_secs*sampling_freq), 'same'));
                    
                end
                
                subplot(no_folders, 2, 2*fo - (2 - ch))
                
                plot(t, beta_blocks_plot(:, :, ch, ty))
                
                axis tight
                
                if ch == 1
                    
                    ylabel(folder)
                    
                    if fo == 1
                        
                        legend(band_labels)
                        
                        title([chan_labels{ch}, ', Proportion High Band Power Per Minute (z-Scored)', long_norms{n}])
                        
                    end
                    
                elseif fo == 1
                    
                    title(chan_labels{ch})
                    
                end
                
                hold on
                
                plot([0; 0], [all_dimensions(@min, beta_blocks_plot(:, :, ch, ty)); all_dimensions(@max, beta_blocks_plot(:, :, ch, ty))], 'r')
                
            end
            
        end
    
        save_as_pdf((n - 1)*2 + ty, [subject_mat(1:(end - length('_subjects.mat'))), '_wt_BP_high', norms{n}, high_type{ty}])
    
    end
    
end

end