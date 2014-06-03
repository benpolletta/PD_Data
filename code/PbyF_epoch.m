function pbf = PbyF_epoch(epoch_name, f_bins)

sampling_freq = 1000;

no_f_bins = length(f_bins) - 1;

P = load(epoch_name);

P_diff = diff(unwrap(P),[],2);

P_diff = angle(exp(sqrt(-1)*P_diff));

F = diff(unwrap(P))/(2*pi*(1/sampling_freq));

F_smooth = nan(size(F));

pbf = nan(1, no_f_bins, 6, 2);

for ch = 1:2
    
    F_smooth(:,ch) = conv(F(:,ch),ones(50,1)/50,'same');
    
    for f = 1:no_f_bins
        
        bin_indicator = f_bins(f) <= F_smooth(:,ch) & F_smooth(:,ch) < f_bins(f+1);
        
        n = sum(bin_indicator); pbf(1, f, 1, ch) = n;
        
        bin_freqs = F_smooth(bin_indicator, ch);
        
        pbf(1, f, 2, ch) = mean(bin_freqs);
        
        pbf(1, f, 3, ch) = std(bin_freqs);
        
        bin_phases = P_diff(bin_indicator);
        
        mrv = sum(exp(sqrt(-1)*bin_phases)); pbf(1, f, 4, ch) = mrv;
        
        pbf(1, f, 5, ch) = abs(mrv)^2 / n;
        
        pbf(1, f, 6, ch) = exp(sqrt(1+4*n+4*(n^2-abs(mrv)^2))-(1+2*n));
        
    end
    
end

save([epoch_name(1:end-4),'_PF.mat'],'P_diff','F_smooth')