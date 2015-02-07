function PD_fix_Spec_pct(subjects_mat, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    suffix = '';
    
else

    suffix = sprintf('%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));

end

load(subjects_mat)

norms = {'', '_norm'}; no_norms = length(norms);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_all_channel_data_dec.mat'], 'sampling_freq')
    
    for n = 1:no_norms
        
        clear All_spec_pct
        
        All_spec = load([subj_name, '_wt', norms{n}, '.mat'], ['Spec', norms{n}]);
        
        All_spec = getfield(All_spec, ['Spec', norms{n}]);
        
        t = (1:size(All_spec, 1))/sampling_freq - basetime;
        
        for ch = 1:2
            
            baseline_mean = mean(abs(All_spec(t <= 0, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
            
            All_spec_pct(:, :, ch) = 100*abs(All_spec(:, :, ch))./(ones(size(All_spec(:, :, ch)))*diag(baseline_mean)) - 100;
            
        end
        
        eval(['Spec', norms{n}, '_pct = All_spec_pct;'])
        
    end
    
    save([subj_name, suffix, '_wt_pct.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_pct', '-v7.3')
    save([subj_name, suffix, '_wt_norm_pct.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_norm_pct', '-v7.3')
        
end