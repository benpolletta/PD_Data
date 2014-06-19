function [MR_vec,p_hist,freqs]=rose_plot(p_vec,f_vec,no_p_bins,f_bins)

% Produces a rose plot, with mean resultant vectors, from time series of
% frequencies and phases.
% INPUTS:
% p_vec = vector of phases, between -pi and pi.
% f_vec = vector of frequencies.
% no_dp_bins = number of bins to divide circle into.
% f_bins = vector of bin edges. Frequencies will be pooled together inside
% each bin.
% OUTPUTS:
% MR_vec = mean resultant vectors for each frequency bin.
% freqs = center frequency for each bin.

p_bins = linspace(-pi,pi,no_p_bins+1)-2*pi/(2*no_p_bins); % Making vector of phase bin edges.
p_bins(1) = pi*(1-1/no_p_bins);   % Taking care of phases close to pi.
phases = (p_bins(1:end-1)+p_bins(2:end))/2; % Making vector of phase bin centers.
phases(1) = -pi; phases(end+1) = -pi; % Taking care of phases close to pi, wrapping vector around end for plotting.

no_f_bins = length(f_bins)-1; % Finding number of frequency bins.

c_order = [linspace(1,0,no_f_bins); abs(linspace(1,0,no_f_bins)-.5); linspace(0,1,no_f_bins)]';

p_hist = zeros(no_p_bins+1,no_f_bins);

freqs = nan(no_f_bins,1);

for f = 1:no_f_bins
    
    freqs(f) = mean(f_bins(f:f+1));
    
    freq_legend{f} = sprintf('%g Hz',freqs(f));
    
    f_phases = p_vec(f_bins(f) <= f_vec & f_vec < f_bins(f+1)); % Finding phases falling within a given frequency bin. 
    
    R = nansum(exp(sqrt(-1)*f_phases));
    
    n = sum(~isnan(f_phases));
    
    MR_pval = exp(sqrt(1+4*n+4*(n^2-abs(R)^2))-(1+2*n));
    
    if MR_pval <= 0.01/no_f_bins && length(f_phases) >= no_p_bins
        
        MR_vec(f) = R/n;
        
        % p_hist(1,f) = sum(p_bins(1) <= f_phases | f_phases < p_bins(2))/length(f_phases); % Normalized (i.e., as a fraction of total phases).
        
        p_hist(1,f) = sum(p_bins(1) <= f_phases | f_phases < p_bins(2)); % Non-normalized.
        
        p_hist(end,f) = p_hist(1,f);
        
        for p = 2:no_p_bins
            
            % p_hist(p,f) = sum(p_bins(p) <= f_phases & f_phases < p_bins(p+1))/length(f_phases); % Normalized (i.e., as a fraction of total phases).
            
            p_hist(p,f) = sum(p_bins(p) <= f_phases & f_phases < p_bins(p+1)); % Non-normalized.
            
        end
        
    else
        
        MR_vec(f) = nan;
        
        p_hist(:,f) = nan;
        
    end
    
end
    
%     rose(phases,no_p_bins)

polar_lim=polar(0,max(max(abs(MR_vec)),max(max(p_hist))),'.'); % Dummy plot to set radial axis length.

set(polar_lim,'Marker','none'); % Dummy plot disappears.

hold on

for f = 1:no_f_bins

    h1(f) = polar(phases',p_hist(:,f),'-');
    
    set(h1(f),'Color',c_order(f,:),'LineWidth',2)
    
    h2(f) = compass(max(max(p_hist))*MR_vec(f)/max(abs(MR_vec)));
    
    set(h2(f),'Color',c_order(f,:),'LineWidth',2)

end

colormap(c_order)

h = colorbar('YTick',(1:no_f_bins)+.5,'YTickLabel',freq_legend);

cbfreeze(h)

% legend(h1,freq_legend)

