function PD_fix_wav_BP(subjects_mat, peak_suffix, freqs, no_cycles, bands, new_bands)

subjects_info = load(subjects_mat);

folders = subjects_info.folders; prefixes = subjects_info.prefixes; basetimes = subjects_info.basetimes;

parfor fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    PD_fix_wav_BP_inner(folder, prefix, peak_suffix, basetime, freqs, no_cycles, bands, new_bands)
        
end

end

function PD_fix_wav_BP_inner(folder, prefix, peak_suffix, basetime, freqs, no_cycles, bands, new_bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, 200);
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    suffix = peak_suffix;
    
else

    suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    suffix = [suffix, peak_suffix];
    
end

if isempty(new_bands)
    
    new_suffix = peak_suffix;
    
else
    
    new_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(new_bands, 1));
    
    new_suffix = [new_suffix, peak_suffix];

end

clear bands

bands = new_bands;

no_bands = size(bands, 1);

norms = {'', '_norm', '_pct'}; no_norms = length(norms); % , '_norm_pct'}; 

subj_name = [folder,'/',prefix];

load([subj_name,'_all_channel_data_dec.mat'])

t = (1:length(PD_dec))/sampling_freq - basetime;

clear BP BP_norm BP_pct

[BP, BP_norm, BP_pct] = deal(nan(length(PD_dec), no_bands, 2));

for n = 1:no_norms
    
    All_spec = load([subj_name, suffix, '_wt', norms{n}, '.mat'], ['Spec', norms{n}]);
    
    All_spec = getfield(All_spec, ['Spec', norms{n}]);
    
    for ch = 1:2
        
        %% Band power.
        
        for b = 1:no_bands
            
            band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
            
            if strcmp(norms{n}, '')
                
                eval(['BP', norms{n}, '(:, b, ch) = sqrt(sum(abs(All_spec(:, band_indices, ch)).^2, 2));'])
                
            else
                
                eval(['BP', norms{n}, '(:, b, ch) = sum(All_spec(:, band_indices, ch), 2);'])
                
            end
            
        end
        
    end
    
    save([subj_name, new_suffix, '_wt_BP.mat'], '-v7.3', 'sampling_freq', 't', 'basetime', 'bands', 'BP', 'BP_norm', 'BP_pct') % , 'BP_norm_pct')
    
end
    
end