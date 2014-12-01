function PD_laser_wav_dec(subjects_mat)

load(subjects_mat)

norms = {'', '_norm'}; no_norms = length(norms);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_all_channel_data_dec.mat'])
    
    load([subj_name, '_wav_laser_artifacts.mat'], 'laser_transitions', 'laser_periods')
    
    laser_transitions = wavelet_spectrogram(laser_transitions, sampling_freq, 1, 3, 0, '');
    
    for n = 1:no_norms
        
        Spec_data = load([subj_name,'_wt.mat'], ['Spec', norms{n}]);
        
        Spec_data = getfield(Spec_data, ['Spec', norms{n}]);
        
        freqs = load([subj_name,'_wt.mat'], 'freqs');
        
        freqs = freqs.freqs;
        
        no_freqs = length(freqs);
        
        Spec_data(laser_transitions > 0) = nan;
        
        pd_blocks = [];
        
        for pd = 1:3
            
            pd_blocks(:, :, pd) = index_to_blocks(laser_periods(:, pd));
            
        end
        
        no_blocks = size(pd_blocks, 1);
        
        Spec_mean = nan(no_blocks, no_freqs, 2, 3);
        
        for pd = 1:3
            
            for b = 1:no_blocks
                
                block_start = max(1, pd_blocks(b, 1, pd));
                
                block_end = min(size(Spec_data, 1), pd_blocks(b, 2, pd));
                
                Spec_mean(b, :, :, pd) = nanmean(Spec_data(block_start:block_end, :, :));
                
            end
            
        end
        
        save([subj_name, '_wav', norms{n}, '_laser_blocks.mat'], 'Spec_mean', 'pd_blocks')
        
    end
    
end
