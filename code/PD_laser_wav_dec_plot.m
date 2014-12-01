function PD_laser_wav_dec_plot(subjects_mat)

load(subjects_mat)

norms = {'', '_norm'}; no_norms = length(norms);

norm_labels = {'Power', 'Power, % Total'};

pd_labels = {'Laser On', 'Post-Laser', 'Pre-Laser'};

for n = 1:no_norms
    
    All_laser_blocks = [];
    
    for fo = 1:length(folders)
        
        folder = folders{fo};
        
        prefix = prefixes{fo};
        
        subj_name = [folder,'/',prefix];
        
        load([subj_name, '_all_channel_data_dec.mat'])
        
        load([subj_name, '_wav_laser_artifacts.mat'], 'laser_transitions', 'laser_periods')
        
        laser_transitions = wavelet_spectrogram(laser_transitions, sampling_freq, 1, 3, 0, '');
        
        freqs = load([subj_name,'_wt.mat'], 'freqs');
        
        freqs = freqs.freqs; plot_indices = freqs >= 0 & freqs <= 50;
        
        load([subj_name, '_wav', norms{n}, '_laser_blocks.mat'])
        
        All_laser_blocks(end + (1:size(Spec_mean, 1)), :, :, :) = Spec_mean;
        
        Sm_size = size(Spec_mean);
        
        pd_mean = reshape(nanmean(Spec_mean), Sm_size(2:end)); pd_mean = permute(pd_mean, [1 3 2]);
        
        pd_se = reshape(nanstd(Spec_mean), Sm_size(2:end))/sqrt(size(Spec_mean, 1)); pd_se = permute(pd_se, [1 3 2]);
        
        figure
        
        for ch = 1:2
            
            subplot(1, 2, ch)
            
            boundedline(freqs(plot_indices), pd_mean(plot_indices, :, ch), prep_for_boundedline(pd_se(plot_indices, :, ch)))
            
            legend(pd_labels)
            
            title([folder, ', ', chan_labels{ch}, ' ', norm_labels{n}])
        
        end
        
    end
    
    Alb_size = size(All_laser_blocks);
    
    all_mean = reshape(nanmean(All_laser_blocks), Alb_size(2:end)); all_mean = permute(all_mean, [1 3 2]);
    
    all_se = reshape(nanstd(All_laser_blocks), Alb_size(2:end))/sqrt(size(All_laser_blocks, 1)); all_se = permute(all_se, [1 3 2]);
    
    figure
    
    for ch = 1:2
        
        subplot(1, 2, ch)
        
        boundedline(freqs(plot_indices), all_mean(plot_indices, :, ch), prep_for_boundedline(all_se(plot_indices, :, ch)))
        
        legend(pd_labels)
        
        title([chan_labels{ch}, ' ', norm_labels{n}])
        
    end
    
end
