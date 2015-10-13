function make_phase_figures(group_prefix, selected_subj_no, c_axis, epoch_secs, f_tol)

close('all')

load([group_prefix, '_subjects.mat'])

no_pds = length(pd_labels); no_chans = length(chan_labels);

time_window = 200; sampling_freq = 500; percent = .5;

no_pds_plotted = 2;

%% Plotting consolidated beta spectrograms.

if ~isempty(epoch_secs)
    
    epoch_label = sprintf('_%ddensest', epoch_secs);
    
else
    
    epoch_label = '';
    
end

figure(2), figure(1)

cons_name = [folders{selected_subj_no}, '/', prefixes{selected_subj_no}, '_2sd_BP_high_',...
    num2str(time_window/sampling_freq), 's_', num2str(percent), 'pct_consolidated'];

for pd = 1:no_pds_plotted

    load([cons_name, epoch_label, '_', pd_labels{pd} '_plot_data.mat'])
    
    for ch = 1:no_chans
        
        if ch == 1
            
            if exist('striatal_id')
            
                chan_id = striatal_id(selected_subj_no);
                
            else
                
                chan_id = find(strcmp(chan_labels, 'Striatum'));
                
            end
            
            chan_title = 'Striatum';

        else
            
            if exist('striatal_id')
            
                chan_id = 3 - striatal_id(selected_subj_no);
                
            else
                
                chan_id = find(~strcmp(chan_labels, 'Striatum'));
                
            end
            
            chan_title = chan_labels{~strcmp(chan_labels, 'Striatum')};
            
        end
            
        subplot(no_pds_plotted, no_chans, (pd - 1)*no_chans + ch) % 2*no_pds_plotted
        
        Spec_hc_chan = Spec_hc_plot(:, :, chan_id);
        
        t_hc_chan = t_hc_plot;
        
        nan_indices = sum(isnan(Spec_hc_chan), 2) == size(Spec_hc_plot, 2);
        
        t_hc_chan(nan_indices) = [];
        
        Spec_hc_chan(nan_indices, :) = [];
        
        if ch == 1
            
            plot_length = min(length(t_hc_chan), 500);
            
        end
        
        plot_start = max(round(length(t_hc_chan)/2) - floor(plot_length/2), 1);
        
        plot_end = min(round(length(t_hc_chan)/2) + floor(plot_length/2) - 1, length(t_hc_chan));
        
        if plot_end > plot_start
            
            imagesc(plot_start:plot_end, 1:50, abs(Spec_hc_chan(plot_start:plot_end, :))')
            
            if pd == 1
                
                title(chan_title, 'FontSize', 20)
                
            end
            
            caxis(c_axis)
            
            axis xy
            
            if ch == 1
                
                ylabel({pd_labels{pd}; 'Freq. (Hz)'}, 'FontSize', 16)
                
            end
            
            if pd == 2
                
                xlabel('Time (s)', 'FontSize', 16)
                
            end
            
            set(gca, 'XTick', [plot_start plot_end], 'XTickLabel', [0 (plot_end - plot_start + 1)*2]/1000, 'FontSize', 16)
            
        end
        
    end
    
end

save_as_eps(1, [group_prefix, '_', folders{selected_subj_no}, '_', num2str(f_tol),'ftol'])

save_as_pdf(1, [group_prefix, '_', folders{selected_subj_no}, '_', num2str(f_tol),'ftol'])

%% Plotting phase angles.

load('bonferroni_count.mat')

no_f_bins = ceil((diff([8 30]) + 1)/f_tol);

f_bins = (mean([8 30]) - f_tol*no_f_bins/2):f_tol:(mean([8 30]) + f_tol*no_f_bins/2);

[MR_mat, conf_mat] = deal(nan(length(f_bins) - 1, no_pds_plotted)); % no_pds));

%figure(1)

[freq_cats, phases] = deal(cell(no_pds_plotted, 1));
    
clear rao_test rayleigh_test zero_test conc_pval angle_pval

for pd = 1:no_pds_plotted % no_pds
    
    clear f_overlap_index mean_f d_phi
    
    [Freqs, Phases] = deal([]);
    
    for ch = 1:no_chans
        
        if ~exist('striatal_id')
            
            if ch == 1
                
                chan_id = find(strcmp(chan_labels, 'Striatum'));
                
            else
                
                chan_id = find(~strcmp(chan_labels, 'Striatum'));
                
            end
            
        else
            
            chan_id = ch;
            
        end
        
        Freq_data = load([group_prefix, '_2sd_8-30Hz_high_', num2str(time_window/sampling_freq), 's_', num2str(percent), 'pct_consolidated_freqs',...
            '_ch', num2str(chan_id), '_', pd_labels{pd}, '.txt']);
    
        if ~isempty(Freq_data)

            Freqs(:, ch) = Freq_data(:, 1);

            Phases(:, ch) = Freq_data(:, 2);

        end
    
    end
    
    f_overlap_index = abs(diff(Freqs, [], 2)) <= f_tol;
    
    mean_f = sum(Freqs, 2)/2;
        
    d_phi = -diff(Phases, [], 2); % Since Striatum will always be first channel, due to chan_id above.
    
    figure
    
    [MR_mat(:, pd), ~, freq_bin_centers, conf_mat(:, pd)] = rose_plot(d_phi(f_overlap_index), mean_f(f_overlap_index), 18, f_bins, bonferroni_count);
    
    freq_cats{pd} = categorize_freq(mean_f(f_overlap_index), f_bins);
    
    phases{pd} = d_phi(f_overlap_index);
    
end

conf_mat = reshape(conf_mat, size(conf_mat, 1), 1, size(conf_mat, 2));

conf_mat = repmat(conf_mat, [1 2 1]);

[~, ~, zero_test, ~, angle_pval] = circular_stats(f_bins, phases, freq_cats, bonferroni_count);

figure(2) % figure(1)

% subplot(2, 1, 2)

h = boundedline(freq_bin_centers', (180/pi)*angle(MR_mat), (180/pi)*conf_mat);

set(h, 'Marker', 's')

axis tight

zero_test(isnan(MR_mat)) = nan;

logical = [(angle_pval < .05) zero_test]; % [(angle_pval < .05/length(angle_pval)) zero_test];

logical(logical == 0) = nan;

add_stars(gca, freq_bin_centers', logical, [1 zeros(1, no_pds_plotted)], [1 0 0; 0 0 1; 0 .5 0])

% axis tight

set(gca, 'NextPlot', 'add', 'ColorOrder', [0 0 1; 0 .5 0])

h_prime = plot(freq_bin_centers', (180/pi)*angle(MR_mat));

set(h_prime, 'Marker', 's')

plot(freq_bin_centers', zeros(length(freq_bin_centers), 1), '--k')

chan_label = ['Striatum - ', chan_labels{~strcmp(chan_labels, 'Striatum')}];

title(['Phase Lag, (', chan_label, ')'], 'FontSize', 20)

ylabel('Degrees', 'FontSize', 16)

xlabel('Frequency (Hz)', 'FontSize', 16)

set(gca, 'FontSize', 16)

legend(h_prime, pd_labels(1:no_pds_plotted))

save_as_eps(2, [group_prefix, '_', num2str(f_tol),'ftol_phase'])

save_as_pdf(2, [group_prefix, '_', num2str(f_tol),'ftol_phase'])

end

%% CATEGORIZE_FREQ
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

%% CIRCULAR_STATS
function [rao_test, rayleigh_test, zero_test, conc_pval, angle_pval] = circular_stats(f_bins, phases, freq_cats, bonferroni_count)

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
rao_test = min(rao_test*bonferroni_count, 1);

rayleigh_test = min(rayleigh_test*bonferroni_count, 1);
        
%% Testing phases against zero.

zero_test = nan(no_f_bins, no_pds);

for f = 1:no_f_bins
    
    for pd = 1:no_pds
        
        phi = phases{pd}(freq_cats{pd} == f); phi(isnan(phi)) = [];
        
        if ~isempty(phi)
            
            zero_test(f, pd) = circ_mtest(phi, 0);
            
        end
        
    end
    
end

% Bonferroni correcting p-values.
zero_test = min(zero_test*bonferroni_count, 1);

%% Testing phases pre- vs. post-infusion.
    
conc_pval = nan(no_f_bins, 1); angle_pval = nan(no_f_bins, 1);

if no_pds > 1
    
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
    conc_pval = min(conc_pval*bonferroni_count, 1); angle_pval = min(angle_pval*bonferroni_count, 1);
    
end

end