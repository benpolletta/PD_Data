function PD_decimate_wav(subjects_mat, epoch_length)

load(subjects_mat)

sampling_freq = 1000;

freqs = 1:200; no_freqs = length(freqs);

bands = [1 4; 4 8; 8 30; 30 100; 100 120; 0 200]; no_bands = size(bands, 1);

norms = {'', '_norm', '_pct'}; no_norms = length(norms);

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
    
    % load([subj_name, '_wt.mat'])
    % 
    % All_spec = nan([size(Spec) 3]);
    % 
    % All_spec(:, :, :, 1) = Spec; All_spec(:, :, :, 2) = Spec_norm; All_spec(:, :, :, 3) = Spec_pct;
    % 
    % clear Spec Spec_norm Spec_pct
    % 
    % All_BP = nan([size(BP) 3]);
    % 
    % All_BP(:, :, :, 1) = BP; All_BP(:, :, :, 2) = BP_norm; All_BP(:, :, :, 3);
    % 
    % clear BP BP_norm BP_pct
        
    t_dec = nan(no_epochs, 1);
    
    Spec_dec = nan(no_epochs, no_freqs, 2, no_norms, 2);
    
    BP_dec = nan(no_epochs, no_bands, 2, no_norms, 2);
    
    for n = 1:no_norms
        
        All_spec = load([subj_name, '_wt.mat'], ['Spec', norms{n}]);
        
        All_BP = load([subj_name, '_wt.mat'], ['BP', norms{n}]);
        
        no_secs = min(basetime + 1500, size(All_spec, 1)/sampling_freq);
        
        t = (1:(no_secs*sampling_freq))/sampling_freq - base_index;
        
        for ch = 1:2
            
            for e = 1:no_epochs
                
                epoch_start = start_index + (e - 1)*epoch_length + 1
                
                epoch_end = start_index + e*epoch_length

		length(t)
                
                t_dec(e) = median(t(epoch_start:epoch_end));
                
                epoch_spec = All_spec(epoch_start:epoch_end, :, ch);
                
                Spec_dec(e, :, ch, n, 1) = median(epoch_spec);
                
                Spec_dec(e, :, ch, n, 2) = mean(epoch_spec);
                
                epoch_BP = All_BP(epoch_start:epoch_end, :, ch);
                
                BP_dec(e, :, ch, n, 1) = median(epoch_BP);
                
                BP_dec(e, :, ch, n, 2) = mean(epoch_BP);
                
            end
            
        end
        
    end
   
    save([subj_name, '_', num2str(epoch_length/sampling_freq), 's_dec_wav.mat'], '-v7.3', 't_dec', 'Spec_dec', 'BP_dec')
    
end
