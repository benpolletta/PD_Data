function PD_beta_blocks_rel_infusion_laser_plots_individual(subject_mat, measure, norm_for_power, no_trials, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

subj_mat_name = subject_mat(1:(end - length('_subjects.mat')));
    
% close('all')

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(chan_labels);

no_pds = length(pd_labels);

pd_indices = cell(no_folders, 1);

if isempty(no_trials)

    trials = cell(no_folders, no_pds);
    
    no_trials = nan(no_folders, no_pds);
    
    for fo = 1:no_folders
        
        folder = folders{fo};
        
        prefix = prefixes{fo};
        
        subj_name = [folder,'/',prefix];
        
        load([subj_name, '_wav_laser_artifacts.mat'])
        
        pd_indices{fo} = logical(laser_periods);
        
        for pd = 1:no_pds
            
            trials{fo, pd} = index_to_blocks(laser_periods(:, pd));
            
            no_trials(fo, pd) = size(trials, 1);
            
        end
        
    end
    
    max_no_trials = all_dimensions(@max, no_trials);
    
else
    
    max_no_trials = no_trials;
    
end

no_secs = max_no_trials*5;
    
load([subj_mat_name, BP_suffix, '_pct_BP_high_laser_', num2str(no_trials), 'trials', measure, norm_for_power, '.mat'])

no_comparisons = nchoosek(no_pds, 2);
    
p_vals = nan(no_comparisons, no_folders, no_bands, no_chans, 2);

p_val_labels = comp_labels(pd_labels);

All_p_vals = nan(no_comparisons, no_bands, no_chans, 2);

format = make_format(5*no_pds + 2*no_comparisons, 'f');

[mean_labels, se_labels, median_labels, q1_labels, q3_labels] = deal({});

for i = 1:no_pds
   
    mean_labels = {mean_labels{:}, 'Mean'};
    
    se_labels = {se_labels{:}, 'S.E.'};
    
    median_labels = {median_labels{:}, 'Median'};
    
    q1_labels = {q1_labels{:}, 'Q1'};
    
    q3_labels = {q3_labels{:}, 'Q3'};
    
end

format = ['%s\t', format];

for b = 1:no_bands
    
    fid = nan(no_chans, 1);
    
    for ch = 1:no_chans
        
        fid(ch) = fopen([subj_mat_name, BP_suffix, '_', num2str(no_trials), '_trials_', short_band_labels{b},...
            '_individual', measure, norm_for_power, '_', chan_labels{ch}, '_stats.txt'], 'w');
        
        fprintf(fid(ch), make_format(1 + 5*no_pds + 2*no_comparisons, 's'), 'Recording', mean_labels{:}, se_labels{:},...
            median_labels{:}, q1_labels{:}, q3_labels{:}, p_val_labels{:});
        
        fprintf(fid(ch), make_format(1 + 5*no_pds + 2*no_comparisons, 's'), 'Recording', pd_labels{:}, pd_labels{:},...
            pd_labels{:}, pd_labels{:}, pd_labels{:}, p_val_labels{:});
        
    end

    All_mean = nan(no_folders, no_pds, no_chans);

    figure
    
    for ch = 1:no_chans
        
        for fo = 1:no_folders
            
            pct_bp_high_for_test = nan(no_secs, no_pds);
            
            for pd = 1:no_pds
                
                pct_bp_high_for_test(1:no_secs, pd) = pct_bp_high(:, fo, pd, b, ch);
                
            end
            
            All_mean(fo, :, ch) = nanmean(pct_bp_high_for_test);
            
            subj_p_vals = run_stats(pct_bp_high_for_test, 'ranksum');
            
            p_vals(:, fo, b, ch, :) = permute(subj_p_vals, [1 3 4 5 2]);
            
            % p_vals = p_vals*all_dimensions(@sum, ~isnan(p_vals(:, b, :, :)));
            
            subplot(no_chans, no_folders, (ch - 1)*no_folders + fo) % subplot(no_bands, 2, (b - 1)*2 + ch)pct_bp_high_for_test = max(pct_bp_high_for_test, eps);
            
            plot_data(pct_bp_high_for_test, fo, strcmp(measure, '_power'), subj_p_vals)
            
            fprintf(fid(ch), format, folders{fo}, nanmean(pct_bp_high_for_test), nanstd(pct_bp_high_for_test)*diag(1./sum(~isnan(pct_bp_high_for_test))),...
                nanmedian(pct_bp_high_for_test), quantile(pct_bp_high_for_test, .25), quantile(pct_bp_high_for_test, .75),...
                reshape(subj_p_vals, 1, 2*no_comparisons));
            
            set(gca, 'XTick', 1:no_pds, 'XTickLabel', pd_labels)
            
            title([folders{fo}, ', ', band_labels{b}])
            
        end
        
        save_as_pdf(gcf, [subj_mat_name, BP_suffix, '_pct_BP_high_laser_', num2str(no_trials), 'trials_',...
            short_band_labels{b}, measure, norm_for_power, '_bar_individual'])
        
    end
    
    %% Stats & figures treating each individual as an observation.
        
    figure
    
    for ch = 1:no_chans
        
        across_p_vals = run_stats(All_mean(:, :, ch), 'ranksum');
        
        All_p_vals(:, b, ch, :) = permute(across_p_vals, [1 3 4 2]);
        
        subplot(1, no_chans, ch)
        
        plot_data(All_mean(:, :, ch), 1, strcmp(measure, '_power'), across_p_vals)
            
        fprintf(fid(ch), format, 'Mean', nanmean(All_mean(:, :, ch)), nanstd(All_mean(:, :, ch))/no_folders,...
            nanmedian(All_mean(:, :, ch)), quantile(All_mean(:, :, ch), .25), quantile(All_mean(:, :, ch), .75),...
            reshape(across_p_vals, 1, 2*no_comparisons));
        
        if no_pds == 2
            
            title({[chan_labels{ch}, ', ', band_labels{b}];...
                ['p-value = ', num2str(across_p_vals(1)), ' (', test_handle, ')']})
            
        else
            
            title([chan_labels{ch}, ', ', band_labels{b}])
            
        end
        
        set(gca, 'XTick', 1:no_pds, 'XTickLabel', pd_labels)
        
    end
    
    save_as_pdf(gcf, [subj_mat_name, BP_suffix, '_pct_BP_high_laser_', num2str(no_trials), 'trials_',...
        short_band_labels{b}, measure, norm_for_power '_bar_individual_avg'])
    
end

end

function pos_bars = get_bar_pos(handle)

    for i = 1:length(handle)

        x = get(get(handle(i), 'children'), 'xdata');

        x = mean(x([1 3],:));

        pos_bars(i,:) = x;

    end

end

function p_vals = run_stats(data_for_test, test_handle)

no_pds = size(data_for_test, 2);

no_comparisons = nchoosek(no_pds, 2);

comparisons = nchoosek(1:no_pds, 2);

p_vals = nan(no_comparisons, 2);

if no_pds == 2 && any(~isnan(data_for_test(:, 1))) && any(~isnan(data_for_test(:, 2)))
    
    if strcmp(test_handle, 't-test')
        
        [~, p_vals(1)] = ttest(data_for_test(:, 1), data_for_test(:, 2), [], 'left'); % 'tail', 'left');
        
    elseif strcmp(test_handle, 'ranksum')
        
        p_vals(1) = ranksum(data_for_test(:, 1), data_for_test(:, 2), 'tail', 'left');
        
    end
    
else
    
    for comp = 1:no_comparisons
        
        if any(~isnan(data_for_test(:, comparisons(comp, 1)))) && any(~isnan(data_for_test(:, comparisons(comp, 2))))
            
            if strcmp(test_handle, 'ranksum')
                
                p_vals(comp, 1) = ranksum(data_for_test(:, comparisons(comp, 1)), data_for_test(:, comparisons(comp, 2)), 'tail', 'left');
                
                p_vals(comp, 2) = ranksum(data_for_test(:, comparisons(comp, 1)), data_for_test(:, comparisons(comp, 2)), 'tail', 'right');
                
            elseif strcmp(test_handle, 't-test')
                
                [~, p_vals(comp, 1)] = ttest(data_for_test(:, comparisons(comp, 1)), data_for_test(:, comparisons(comp, 2)), [], 'left');
                
                [~, p_vals(comp, 2)] = ttest(data_for_test(:, comparisons(comp, 1)), data_for_test(:, comparisons(comp, 2)), [], 'right');
                
            end
            
        end
        
    end
    
end

end

function p_val_labels = comp_labels(pd_labels)

no_pds = length(pd_labels);

no_comparisons = nchoosek(no_pds, 2);

comparisons = nchoosek(1:no_pds, 2);

p_val_labels = cell(2*no_comparisons);

for comp = 1:no_comparisons
    
    p_val_labels{comp} = [pd_labels{comparisons(comp, 1)}, '<', pd_labels{comparisons(comp, 2)}];
    
    p_val_labels{no_comparisons + comp} = [pd_labels{comparisons(comp, 1)}, '>', pd_labels{comparisons(comp, 2)}];
    
end

end

function plot_data(data, fo, power_flag, p_vals)

no_pds = size(data, 2);

if power_flag
    
    mn = nanmean(data);
    
    sd = nanstd(data)*diag(1./sqrt(sum(~isnan(data))));
    
    if length(mn) == 1 % Stupid Matlab can't tell the difference between x & y and y & width.
        
        mn = [mn nan]; sd = [sd nan];
        
    end
    
    h = barwitherr(sd, mn, 0.6);
    
    if fo == 1
        
        ylabel('Beta Power')
        
    end
    
else
    
    data = max(data, eps);
    
    mn = nanmean(log(data));
    
    sd = nanstd(log(data))*diag(1./sqrt(sum(~isnan(log(data)))));
    
    if length(mn) == 1 % Stupid Matlab can't tell the difference between x & y and y & width.
        
        mn = [mn mn - 5]; sd = [sd nan];
        
    end
    
    h = barwitherr(sd, mn, 0.6, 'BaseValue', min(nanmean(log(data))) - 5);
    
    if fo == 1
        
        ylabel('Beta Density')
    
    end
        
end

box off
            
xlim([0 (no_pds + 1)])

no_comparisons = nchoosek(no_pds, 2);

comparisons = nchoosek(1:no_pds, 2);

bar_pos = get_bar_pos(h);

for comp_type = 1:2
    
    bar_pairs = {};
    
    for comp = 1:no_comparisons
        
        if p_vals(comp, comp_type) < .05
            
            bar_pairs = {bar_pairs{:}, [bar_pos(comparisons(comp, 1)), bar_pos(comparisons(comp, 2))]};
            
        end
        
    end
    
    sigstar(bar_pairs, p_vals(p_vals(:, comp_type) < .05, comp_type))
    
end

end

%% Boxplots.
% 
% [p_greater, p_less] = deal(nan(no_bands, no_comparisons, no_chans));
% 
% for b = 1:no_bands
%     
%     figure
%     
%     for ch = 1:no_chans
%         
%         pct_bp_high_for_test = nan(no_folders*no_secs, no_pds);
%         
%         for pd = 1:no_pds
%             
%             pct_bp_high_for_test(1:no_folders*no_secs, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*no_secs, 1);
%             
%         end
%         
%         for comp = 1:no_comparisons
%         
%             p_greater(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'left');
%         
%             p_less(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'right');
%             
%         end
%         
%         subplot(1, no_chans, ch), % subplot(no_bands, 2, (b - 1)*2 + ch)
%         
%         boxplot(pct_bp_high_for_test, 'labels', pd_labels)
%         
%         ylabel('High Beta Density')
%         
%         % if b == 1
%             
%         title([chan_labels{ch}, ', Boxplot of High ', band_labels{b}, ' Density (Per Trial)'])
%         
%         % else
%         % 
%         %     title([band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
%         % 
%         % end
%         
%     end
% 
%     save_as_pdf(gcf, [subj_mat_name, '_pct_BP_high_laser_trials_', short_band_labels{b}, '_boxplot'])
%     
% end

%% Histograms.
% 
% [p_greater, p_less] = deal(nan(no_bands, no_comparisons, no_chans));
% 
% for b = 1:no_bands
%     
%     figure
%     
%     for ch = 1:no_chans
%         
%         pct_bp_high_for_test = nan(no_folders*no_secs, no_pds);
%         
%         for pd = 1:no_pds
%             
%             pct_bp_high_for_test(1:no_folders*no_secs, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*no_secs, 1);
%             
%         end
%         
%         for comp = 1:no_comparisons
%         
%             p_greater(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'left');
%         
%             p_less(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'right');
%             
%         end
%         
%         subplot(1, no_chans, ch) % subplot(no_bands, 2, (b - 1)*2 + ch)
%         
%         range = [all_dimensions(@min, pct_bp_high_for_test) all_dimensions(@max, pct_bp_high_for_test)];
%         
%         range = max(eps, range);
%         
%         [bin_edges, bin_centers] = make_bins(range(1), range(2), 50, 'log'); % 10, ''); %
%         
%         for pd = 1:no_pds
%             
%             [h, ~] = histc(pct_bp_high_for_test(:, pd), bin_edges);
%             
%             % h = max(eps, h);
%             
%             loglog(bin_centers, h(1:(end - 1))/sum(h(1:(end - 1))), pd_colors{pd})
%             
%             hold on
%             
%         end
%         
%         legend(pd_labels)
%         
%         xlabel('High Beta Density')
%         
%         ylabel('Proportion Observed')
%             
%         title([chan_labels{ch}, ', Histogram of High ', band_labels{b}, ' Density (Per Trial)'])
%         
%         % if b == 1
%         % 
%         %     title([chan_labels{ch}, ', ', band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
%         % 
%         % else
%         % 
%         %     title([band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
%         % 
%         % end
%         
%     end
% 
%     save_as_pdf(gcf, [subj_mat_name, '_pct_BP_high_laser_trials_', short_band_labels{b}, '_hist'])
%     
% end