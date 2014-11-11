function PD_pmtm(subjects_mat, epoch_length)

load(subjects_mat)

sampling_freq = 1000;

bands = [1 4; 4 8; 8 30; 30 100; 100 120; 0 200]; no_bands = size(bands, 1);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    base_index = basetime*sampling_freq;
    
    no_pre_epochs = floor(base_index/epoch_length);
    
    start_index = base_index - no_pre_epochs*epoch_length;
    
    no_post_epochs = floor(1500*sampling_freq/epoch_length);
    
    no_epochs = no_pre_epochs + no_post_epochs;
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    test_data = PD_dec(1:epoch_length, 1);
    
    [~, freqs] = pmtm(test_data, 1*epoch_length/sampling_freq, [], sampling_freq);
    
    no_freqs = length(freqs);
    
    clear Spec Spec_norm Spec_pct BP BP_norm BP_pct
    
    epoch_no = nan(no_epochs, 1);
    
    [Spec, Spec_norm, Spec_pct] = deal(nan(no_epochs, no_freqs, 2));
        
    [BP, BP_norm, BP_pct] = deal(nan(no_epochs, no_bands, 2));
        
    for e = 1:no_epochs
        
        epoch_no(e) = e - no_pre_epochs - 1;
        
        epoch_start = start_index + (e - 1)*epoch_length + 1;
        
        epoch_end = start_index + e*epoch_length;
        
        epoch_data = PD_dec(epoch_start:epoch_end, :);
        
        for ch = 1:2
            
            data = epoch_data(:, ch);
            
            data_hat = pmtm(data, 1*epoch_length/sampling_freq);
            
            Spec(e, :, ch) = data_hat;
            
            %% Band power.
            
            for b = 1:no_bands
                
                band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
                
                BP(e, b, ch) = sum(abs(data_hat(band_indices)), 2);
                
            end
            
        end
        
    end
    
    for ch = 1:2
        
        %% Normalize by total power.
        
        BP_norm(:, :, ch) = BP(:, :, ch)./repmat(BP(:, end, ch), 1, no_bands);
        
        Spec_norm(:, :, ch) = Spec(:, :, ch)./repmat(BP(:, end, ch), 1, no_freqs);
        
        %% Baseline normalize.
        
        baseline_mean = mean(abs(Spec_norm(t <= basetime, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;
        
        baseline_BP = mean(BP_norm(t <= basetime, :, ch));
        
        BP_pct(:, :, ch) = 100*BP_norm(:, :, ch)./ones(size(BP_norm(:, :, ch)))*diag(baseline_BP) - 100;
        
    end
        
    save([subj_name, '_', num2str(epoch_length/sampling_freq),'s_epoch_pmtm.mat'], '-v7.3', 'epoch_no', 'Spec', 'Spec_norm', 'Spec_pct', 't', 'basetime', 'freqs', 'BP', 'BP_norm', 'BP_pct')
        
end