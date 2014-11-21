function PD_pmtm(subjects_mat, epoch_length)

load(subjects_mat)

sampling_freq = 1000;

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 500]; no_bands = size(bands, 1);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    basetime = basetimes(fo);
    
    base_index = basetime*sampling_freq;
    
    no_pre_epochs = floor(base_index/epoch_length);
    
    start_index = base_index - no_pre_epochs*epoch_length;
    
    no_post_epochs = floor((size(PD_dec, 1) - basetime*sampling_freq)/epoch_length); %floor(min(1500, size(PD_dec, 1)/sampling_freq - basetime)*sampling_freq/epoch_length);
    
    no_epochs = no_pre_epochs + no_post_epochs;
    
    test_data = PD_dec(1:epoch_length, 1);
    
    [~, freqs] = pmtm(test_data, 1*epoch_length/sampling_freq, [], sampling_freq);
    
    no_freqs = length(freqs);
    
    clear Spec Spec_norm Spec_pct Spec_norm_pct BP BP_norm BP_pct BP_norm_pct
    
    [epoch_no, t] = deal(nan(no_epochs, 1));
    
    [Spec, Spec_norm, Spec_pct, Spec_norm_pct] = deal(nan(no_epochs, no_freqs, 2));
        
    [BP, BP_norm, BP_pct, BP_norm_pct] = deal(nan(no_epochs, no_bands, 2));
        
    for e = 1:no_epochs
        
        epoch_no(e) = e - no_pre_epochs - 1;
        
        epoch_start = start_index + (e - 1)*epoch_length + 1;
        
        epoch_end = start_index + e*epoch_length;
        
        t(e) = mean([epoch_start epoch_end])/sampling_freq - basetime;
        
        epoch_data = PD_dec(epoch_start:epoch_end, :);
        
        for ch = 1:2
            
            data = epoch_data(:, ch);
            
            data_hat = pmtm(data, 1*epoch_length/sampling_freq);
            
            Spec(e, :, ch) = data_hat;
            
            %% Band power.
            
            for b = 1:no_bands
                
                band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
                
                BP(e, b, ch) = sum(abs(data_hat(band_indices)));
                
            end
            
        end
        
    end
    
    for ch = 1:2
        
        %% Baseline normalize.
        
        baseline_mean = mean(abs(Spec(t < 0, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_pct(:, :, ch) = 100*abs(Spec(:, :, ch))./(ones(size(Spec(:, :, ch)))*diag(baseline_mean)) - 100;
        
        % baseline_BP = mean(BP(t < 0, :, ch));
        % 
        % BP_pct(:, :, ch) = 100*BP(:, :, ch)./(ones(size(BP(:, :, ch)))*diag(baseline_BP)) - 100;
        
        %% Normalize by total power.
        
        Spec_norm(:, :, ch) = Spec(:, :, ch)./repmat(BP(:, end, ch), 1, no_freqs);
        
        % BP_norm(:, :, ch) = BP(:, :, ch)./repmat(BP(:, end, ch), 1, no_bands);
        
        %% Baseline normalize percent of total power.
        
        baseline_mean = mean(abs(Spec_norm(t < 0, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_norm_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;
         
        % baseline_BP = mean(BP_norm(t < 0, :, ch));
        % 
        % BP_norm_pct(:, :, ch) = 100*BP_norm(:, :, ch)./(ones(size(BP_norm(:, :, ch)))*diag(baseline_BP)) - 100;
        
        %% Band power.
        
        for b = 1:no_bands
            
            band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
            
            BP_pct(:, b, ch) = sum(Spec_pct(:, band_indices, ch), 2);
            
            BP_norm(:, b, ch) = sum(Spec_norm(:, band_indices, ch), 2);
            
            BP_norm_pct(:, b, ch) = sum(Spec_norm_pct(:, band_indices, ch), 2);
            
        end
        
    end
        
    save([subj_name, '_', num2str(epoch_length/sampling_freq),'s_epoch_pmtm.mat'], '-v7.3', 'epoch_no', 'Spec', 'Spec_norm', 'Spec_pct', 'Spec_norm_pct', 't', 'basetime', 'freqs', 'BP', 'BP_norm', 'BP_pct', 'BP_norm_pct')
        
end