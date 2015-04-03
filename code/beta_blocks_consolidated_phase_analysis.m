function beta_blocks_consolidated_phase_analysis(subject_mat, time_window, percent, norm, band_index, f_tol, freqs, no_cycles, bands)

PD_struct = PD_initialize(subject_mat);

if strcmp(PD_struct.chan_labels{2}, 'Striatum'), 
    
    dphi_multiplier = 1; 
    
    chan_label = [PD_struct.chan_labels{2}, ' - ', PD_struct.chan_labels{1}];

else
    
    dphi_multiplier = -1;
    
    chan_label = [PD_struct.chan_labels{1}, ' - ', PD_struct.chan_labels{2}];

end

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_f_bins = ceil((diff(bands(band_index, :)) + 1)/f_tol);

f_bins = (mean(bands(band_index, :)) - f_tol*no_f_bins/2):f_tol:(mean(bands(band_index, :)) + f_tol*no_f_bins/2);

[r, c] = subplot_size(PD_struct.no_folders);

for fo = 1:PD_struct.no_folders
    
    folder = PD_struct.folders{fo};
    
    prefix = PD_struct.prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    cons_name = [subj_name, BP_suffix, '_2sd_BP_high_',...
        num2str(time_window/PD_struct.sampling_freq), 's_', num2str(percent), 'pct_consolidated'];
    
    % cons_title = [folder, ', High Power, ', num2str(percent), ' Percent of ', num2str(time_window/PD_struct.sampling_freq), ' s'];
    
    load([cons_name, '.mat'], 'BP_high_consolidated')
    
    load([subj_name, '_all_channel_data_dec.mat'])
    
    figure(1)
    
    no_phase_bins = 18;
    
    [MR_mat, conf_mat] = deal(nan(length(f_bins) - 1, PD_struct.no_pds));
    
    [freq_cats, phases] = deal(cell(PD_struct.no_pds, 1));
    
    clear mean_f d_phi rao_test zero_test rayleigh_test conc_pval angle_pval
    
    for pd = 1:PD_struct.no_pds    
    
        load([cons_name, '_', short_band_labels{band_index}, '_', PD_struct.pd_labels{pd}, '_freqs', norm, '.mat'])
        
        f_overlap_index = abs(diff(Freqs_high_beta, [], 2)) <= f_tol;
        
        mean_f = sum(Freqs_high_beta, 2)/2;
        
        d_phi = dphi_multiplier*diff(Phases_high_beta, [], 2);
        
        figure
        
        [MR_mat(:, pd), ~, freq_bin_centers, conf_mat(:, pd)] = rose_plot(d_phi(f_overlap_index), mean_f(f_overlap_index), no_phase_bins, f_bins);
        
        freq_cats{pd} = categorize_freq(mean_f(f_overlap_index), f_bins);
        
        phases{pd} = d_phi(f_overlap_index);
        
    end

    [~, ~, zero_test, ~, angle_pval] = circular_stats(f_bins, phases, freq_cats);
    
    conf_mat = reshape(conf_mat, size(conf_mat, 1), 1, size(conf_mat, 2));
    
    conf_mat = repmat(conf_mat, [1 2 1]);
    
    figure(1)
    
    subplot(r, c, fo)
    
    h = boundedline(freq_bin_centers', (180/pi)*angle(MR_mat), (180/pi)*conf_mat);
    
    set(h, 'Marker', 's')
    
    zero_test(isnan(MR_mat)) = nan;
    
    logical = [(angle_pval < .05/length(angle_pval)) zero_test];
    
    logical(logical == 0) = nan;

    add_stars(gca, freq_bin_centers', logical, [1 0 0], [1 0 0; 0 0 1; 0 .5 0])
   
    axis tight
    
    hold on
    
    plot(freq_bin_centers', zeros(length(freq_bin_centers), 1), '--k')
    
    title(['Phase Lag (', chan_label, ')'])
    
    if mod(fo, c) == 1
    
        ylabel('Degrees')
    
        if floor(fo/r) == 0
           
            legend(PD_struct.pd_labels)
            
        end
        
    end
    
    if floor(fo/r) == c
        
        xlabel('Frequency (Hz)')
        
    end
    
end

save_as_pdf(1, [PD_struct.subj_prefix, BP_suffix, '_2sd_BP_high_', short_band_labels{band_index}, '_',...
        num2str(time_window/PD_struct.sampling_freq), 's_', num2str(percent), 'pct_consolidated_', num2str(f_tol), 'ftol_phase_individual'])
    
close('all')
    
[MR_mat, conf_mat] = deal(nan(length(f_bins) - 1, PD_struct.no_pds));

figure(1)

[freq_cats, phases] = deal(cell(PD_struct.no_pds, 1));
    
clear rao_test rayleigh_test zero_test conc_pval angle_pval

for pd = 1:2
    
    clear f_overlap_index mean_f d_phi
    
    [Freqs, Phases] = deal([]);
    
    for ch = 1:2
        
        Freq_data = load([PD_struct.subj_prefix, BP_suffix, '_2sd_', short_band_labels{band_index}, '_high_',...
            num2str(time_window/PD_struct.sampling_freq), 's_', num2str(percent), 'pct_consolidated_freqs',...
            norm, '_ch', num2str(ch), '_', PD_struct.pd_labels{pd}, '.txt']);
    
        if ~isempty(Freq_data)

            Freqs(:, ch) = Freq_data(:, 1);

            Phases(:, ch) = Freq_data(:, 2);

        end
    
    end
    
    f_overlap_index = abs(diff(Freqs, [], 2)) <= f_tol;
    
    mean_f = sum(Freqs, 2)/2;
    
    d_phi = dphi_multiplier*diff(Phases, [], 2);
    
    figure
    
    [MR_mat(:, pd), ~, freq_bin_centers, conf_mat(:, pd)] = rose_plot(d_phi(f_overlap_index), mean_f(f_overlap_index), no_phase_bins, f_bins);
    
    freq_cats{pd} = categorize_freq(mean_f(f_overlap_index), f_bins);
    
    phases{pd} = d_phi(f_overlap_index);
    
end

[~, ~, zero_test, ~, angle_pval] = circular_stats(f_bins, phases, freq_cats);

figure(1)

h = boundedline(freq_bin_centers', (180/pi)*angle(MR_mat), (180/pi)*conf_mat);

set(h, 'Marker', 's')

zero_test(isnan(MR_mat)) = nan;

logical = [(angle_pval < .05/length(angle_pval)) zero_test];

logical(logical == 0) = nan;

add_stars(gca, freq_bin_centers', logical, [1 0 0], [1 0 0; 0 0 1; 0 .5 0])

axis tight

hold on

plot(freq_bin_centers', zeros(length(freq_bin_centers), 1), '--k')

title(['Phase Lag (', chan_label, ')'])

ylabel('Degrees')

xlabel('Frequency (Hz)')

legend(PD_struct.pd_labels)

save_as_pdf(1, [PD_struct.subj_prefix, BP_suffix, '_2sd_BP_high_', short_band_labels{band_index}, '_',...
        num2str(time_window/PD_struct.sampling_freq), 's_', num2str(percent), 'pct_consolidated_', num2str(f_tol), 'ftol_phase_individual_avg'])
    
end

function F_c = categorize_freq(F, f_bins)
    
    [r, c] = size(F);
    
    F_c = zeros(r, c);
    
    no_f_bins = length(f_bins) - 1;
    
    for col = 1:c
        
        for f = 1:no_f_bins
            
            F_bin = F(:, col) >= f_bins(f) & F(:, col) < f_bins(f + 1);
            
            F_c(:, col) = F_c(:, col) + f*F_bin;
            
        end
        
    end
    
end

function [rao_test, rayleigh_test, zero_test, conc_pval, angle_pval] = circular_stats(f_bins, phases, freq_cats)

no_pds = length(phases);

no_f_bins = length(f_bins) - 1;

%% Testing phases for uniformity.

rao_test = nan(no_f_bins, 2);

rayleigh_test = nan(no_f_bins, 2);

for f = 1:no_f_bins
    
    for pd = 1:no_pds
        
        phi = phases{pd}(freq_cats{pd} == f); phi(isnan(phi)) = [];
        
        if ~isempty(phi)
            
            rao_test(f, pd) = circ_raotest(phi);
            
            rayleigh_test(f, pd) = circ_rtest(phi);
            
        end
        
    end
    
end

% Bonferroni correcting p-values.
rao_test = min(rao_test*2*no_f_bins, 1);

rayleigh_test = min(rayleigh_test*2*no_f_bins, 1);
        
%% Testing phases against zero.

zero_test = nan(no_f_bins, 2);

for f = 1:no_f_bins
    
    for pd = 1:2
        
        phi = phases{pd}(freq_cats{pd} == f); phi(isnan(phi)) = [];
        
        if ~isempty(phi)
            
            zero_test(f, pd) = circ_mtest(phi, 0);
            
        end
        
    end
    
end

% Bonferroni correcting p-values.
zero_test = min(zero_test*2*no_f_bins, 1);

%% Testing phases pre- vs. post-infusion.

if no_pds > 1
    
    conc_pval = nan(no_f_bins, 1); angle_pval = nan(no_f_bins, 1);
    
    for f = 1:no_f_bins
            
        if any(~isnan(phases{1})) && any(~isnan(phases{2}))
        
            phi_pre = phases{1}(freq_cats{1} == f); phi_pre(isnan(phi_pre)) = [];
            
            phi_post = phases{2}(freq_cats{2} == f); phi_post(isnan(phi_post)) = [];
            
            conc_pval(f) = circ_ktest(phi_pre, phi_post);
            
            if size(phi_pre, 1) > 0 && size(phi_post, 1) > 0
                
                angle_pval(f) = circ_wwtest(phi_pre, phi_post);
                
            else
                
                angle_pval(f) = 1;
                
            end
            
        else
            
            conc_pval(f) = 1;
            
            angle_pval(f) = 1;
            
        end
        
    end
    
    % Bonferroni correcting p-values.
    conc_pval = min(conc_pval*no_f_bins, 1); angle_pval = min(angle_pval*no_f_bins, 1);
    
end

end