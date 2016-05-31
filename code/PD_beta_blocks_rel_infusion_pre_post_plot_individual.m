function PD_beta_blocks_rel_infusion_pre_post_plot_individual(subject_mat, peak_suffix, epoch_secs, pd_handle, test_handle, freqs, no_cycles, bands, folder_left_out)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    % freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = ['', peak_suffix];
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end

if ~isempty(folder_left_out)
    
    folder_flag = ['_no_', folder_left_out];
   
else
    
    folder_flag = '';
    
end
    
close('all')

load(subject_mat)

no_folders = length(folders);

% if exist([folders{1}, '/', prefixes{1}, '_wt.mat'])
%     
%     load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')
%     
% else
%     
%     sampling_freq = 500;
%     
% end

no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(chan_labels);

no_pds = length(pd_labels);

% pd_colors = {'g', 'r'};

% norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);
% 
% long_norms = {'', ', Increase Over Baseline Power', ', % Total Power', ', Increase in % Total Power Over Baseline'};

% high_type = {'', '_cum'}; no_types = length(high_type);

subj_mat_name = subject_mat(1:(end - length('_subjects.mat')));
    
load([subj_mat_name, BP_suffix '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, '.mat'])

power_flag = exist('BP_sec', 'var') || strcmp(pd_handle, '_power') || strcmp(pd_handle, '_pct_power');

if exist('BP_sec', 'var')
    
    pct_bp_high = BP_sec;
    
end

if no_pds > 1
    
    no_comparisons = nchoosek(no_pds, 2);
    
    p_val_labels = comp_labels(pd_labels);
    
else
    
    no_comparisons = 0;
    
    p_val_labels = {};
    
end

p_vals = nan(no_comparisons, no_folders, no_bands, no_chans, 2);
    
All_p_vals = nan(no_comparisons, no_bands, no_chans, 2);

format = make_format(5*no_pds + 2*no_comparisons, 'f');
    
format = ['%s\t', format];

stat_labels = {'Recording'};

stats = {'Mean', 'S.E.', 'Median', 'Q1', 'Q3'};

for s = 1:length(stats)
    
    for p = 1:no_pds
        
        stat_labels = {stat_labels{:}, stats{s}};
        
    end
    
end

%% Individual bar plots.

for b = 1:no_bands
    
    fid = nan(no_chans, 1);
    
    for ch = 1:no_chans
        
        fid(ch) = fopen([subj_mat_name, BP_suffix, '_', num2str(epoch_secs/60), '_mins_', short_band_labels{b},...
            '_individual_', test_handle, pd_handle, '_', chan_labels{ch}, '_stats.txt'], 'w');
        
        fprintf(fid(ch), make_format(1 + 5*no_pds + 2*no_comparisons, 's'), stat_labels{:}, p_val_labels{:}); %'Recording', 'Mean', 'Mean', 'S.E.', 'S.E.',...
            % 'Median', 'Median', 'Q1', 'Q1', 'Q3', 'Q3', p_val_labels{:});
        
        fprintf(fid(ch), make_format(1 + 5*no_pds, 's'), '', pd_labels{:}, pd_labels{:},...
            pd_labels{:}, pd_labels{:}, pd_labels{:});
        
    end
    
    figure
    
    All_mean = nan(no_folders, no_pds, no_chans);
    
    for fo = 1:no_folders
        
        folder = folders{fo};
        
        prefix = prefixes{fo};
        
        subj_name = [folder, '/', prefix];
        
        for ch = 1:no_chans
            
            pct_bp_high_for_test = nan(epoch_secs, no_pds);
            
            for pd = 1:no_pds
                
                pct_bp_high_for_test(:, pd) = pct_bp_high(:, fo, pd, b, ch);
                
            end
            
            if no_comparisons > 0
                
                subj_p_vals = run_stats(pct_bp_high_for_test, test_handle);
            
                p_vals(:, fo, b, ch, :) = permute(subj_p_vals, [1 3 4 5 2]);
                
            else
               
                subj_p_vals = [];
                
            end
            
            % p_val = p_val*all_dimensions(@sum, ~isnan(p_val))*no_folders;
            
            subplot(no_chans, no_folders, (ch - 1)*no_folders + fo), % subplot(no_bands, 2, (b - 1)*2 + ch)
            
            plot_data(pct_bp_high_for_test, fo, power_flag, subj_p_vals, [chan_labels{ch}, ', ', short_band_labels{b}])
            
            fprintf(fid(ch), format, folder, nanmean(pct_bp_high_for_test), nanstd(pct_bp_high_for_test)/sqrt(epoch_secs),...
                nanmedian(pct_bp_high_for_test), quantile(pct_bp_high_for_test, .25), quantile(pct_bp_high_for_test, .75),...
                reshape(subj_p_vals, 1, 2*no_comparisons));
            
            if ~strcmp(folders{fo}, folder_left_out)
                
                All_mean(fo, :, ch) = nanmean(pct_bp_high_for_test);
                
            end
            
            if no_pds == 2
                
                if fo == 1
                    
                    title({[num2str(epoch_secs/60), ' Mins Densest High Power'];...
                        [folder, ', ', test_handle];['p=', num2str(p_vals(:, fo, b, ch))]})
                    
                else
                    
                    title({[folder, ','];['p=', num2str(p_vals(:, fo, b, ch))]})
                    
                end
                
            else
                
                title(folder)
                
            end
            
            set(gca, 'XTick', 1:no_pds, 'XTickLabel', pd_labels)
            
        end
    
        save([subj_name, BP_suffix, '_pct_BP_high_pvals'], 'p_vals')
        
    end
        
    save_as_pdf(gcf, [subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{b},...
        '_barplot_individual_', test_handle, pd_handle])
    
    %% Stats & figures treating each individual as an observation.
        
    figure
    
    for ch = 1:no_chans
        
        if no_comparisons > 0
            
            across_p_vals = run_stats(All_mean(:, :, ch), [test_handle, '_avg']);
            
            All_p_vals(:, b, ch, :) = permute(across_p_vals, [1 3 4 2]);
            
        else
            
            across_p_vals = [];
            
        end
        
        subplot(1, no_chans, ch)
        
        plot_data(All_mean(:, :, ch), 1, power_flag, across_p_vals, [chan_labels{ch}, ', ', short_band_labels{b}])
            
        fprintf(fid(ch), format, 'Mean', nanmean(All_mean(:, :, ch)), nanstd(All_mean(:, :, ch))/sqrt(no_folders),...
            nanmedian(All_mean(:, :, ch)), quantile(All_mean(:, :, ch), .25), quantile(All_mean(:, :, ch), .75),...
            reshape(across_p_vals, 1, 2*no_comparisons));
        
        if no_pds == 2
            
            title({[chan_labels{ch}, ', ', num2str(epoch_secs/60), ' Mins Densest High Power, ', test_handle];...
                ['p = ', num2str(across_p_vals(1))]})
            
        else
            
            title([chan_labels{ch}, ', ', num2str(epoch_secs/60), ' Mins Densest High Power'])
            
        end
        
        set(gca, 'XTick', 1:no_pds, 'XTickLabel', pd_labels)
        
    end
    
    save_as_pdf(gcf, [subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{b},...
        '_barplot_individual_avg', folder_flag, '_', test_handle, pd_handle])

end

end

function pos_bars = get_bar_pos(handle)

    for i = 1:length(handle)

        x = get(get(handle(i), 'children'), 'xdata');

        x = mean(x([1 3],:));

        pos_bars(i,:) = x;

    end

end

function [p_vals] = run_stats(data_for_test, test_handle)

no_pds = size(data_for_test, 2);

no_comparisons = nchoosek(no_pds, 2);

comparisons = nchoosek(1:no_pds, 2);

p_vals = nan(no_comparisons, 2);

if no_pds == 2 && any(~isnan(data_for_test(:, 1))) && any(~isnan(data_for_test(:, 2)))
    
    if strcmp(test_handle, 't-test')
        
        [~, p_vals(1)] = ttest2(data_for_test(:, 1), data_for_test(:, 2), 'tail', 'left'); % 'tail', 'left');
        
    elseif strcmp(test_handle, 't-test_avg')
        
        [~, p_vals(1)] = ttest(data_for_test(:, 1), data_for_test(:, 2), [], 'left'); % 'tail', 'left');
        
    elseif strcmp(test_handle, 'ranksum') || strcmp(test_handle, 'ranksum_avg')
        
        p_vals(1) = ranksum(data_for_test(:, 1), data_for_test(:, 2), 'tail', 'left');
        
    end
    
else
    
    for comp = 1:no_comparisons
        
        if any(~isnan(data_for_test(:, comparisons(comp, 1)))) && any(~isnan(data_for_test(:, comparisons(comp, 2))))
            
            if strcmp(test_handle, 'ranksum') || strcmp(test_handle, 'ranksum_avg')
                
                p_vals(comp, 1) = ranksum(data_for_test(:, comparisons(comp, 1)), data_for_test(:, comparisons(comp, 2)), 'tail', 'left');
                
                p_vals(comp, 2) = ranksum(data_for_test(:, comparisons(comp, 1)), data_for_test(:, comparisons(comp, 2)), 'tail', 'right');
                
            elseif strcmp(test_handle, 't-test')
                
                [~, p_vals(comp, 1)] = ttest2(data_for_test(:, comparisons(comp, 1)), data_for_test(:, comparisons(comp, 2)), 'tail', 'left');
                
                [~, p_vals(comp, 2)] = ttest2(data_for_test(:, comparisons(comp, 1)), data_for_test(:, comparisons(comp, 2)), 'tail', 'right');
                
            elseif strcmp(test_handle, 't-test_avg')
                
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

function plot_data(data, fo, power_flag, p_vals, y_label)

no_pds = size(data, 2);

if power_flag
    
    mn = nanmean(data);
    
    sd = nanstd(data)*diag(1./sqrt(sum(~isnan(data))));
    
    if length(mn) == 1 % Stupid Matlab can't tell the difference between x & y and y & width.
        
        mn = [mn nan]; sd = [sd nan];
        
    end
    
    h = barwitherr(sd, mn, 0.6);
    
    if fo == 1
        
        ylabel({y_label; 'Power'})
            
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
        
        ylabel({y_label; 'Density'})
    
    end
        
end

box off
            
xlim([0 (no_pds + 1)])

if no_pds > 1
    
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

end

%% Group boxplots.
%
% p_val = nan(no_bands, no_chans);
% 
% for b = 1:no_bands
%     
%     figure
%     
%     for ch = 1:no_chans
%         
%         pct_bp_high_for_test = nan(no_folders*epoch_secs, no_pds);
%         
%         for pd = 1:no_pds
%             
%             pct_bp_high_for_test(:, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*epoch_secs, 1);
%             
%         end
%         
%         if no_pds == 2 %for comp = 1:no_comparisons
%         
%             p_val(b, ch) = ranksum(pct_bp_high_for_test(:, 1), pct_bp_high_for_test(:, 2), 'tail', 'left');
%             
%         end
%         
%         subplot(1, no_chans, ch), % subplot(no_bands, 2, (b - 1)*2 + ch)
%         
%         boxplot(pct_bp_high_for_test, 'labels', pd_labels)
%         
%         ylabel('High Beta Density')
%         
%         % if b ==
%             
%             title({[chan_labels{ch}, ','];[num2str(epoch_secs/60), ' Minutes of Densest High Power'];...
%                 [band_labels{b}, ', p-value = ', num2str(p_val(b, ch)), ' (ranksum test)']})
%         
%         % else
%         % 
%         %     title([band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
%         % 
%         % end
%         
%     end
% 
%     save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{b}, '_boxplot'])
%     
% end

%% Group histograms.
% 
% p_val = nan(no_bands, no_chans);
% 
% for b = 1:no_bands
%     
%     figure
%     
%     for ch = 1:no_chans
%         
%         pct_bp_high_for_test = nan(no_folders*epoch_secs, 2);
%         
%         for pd = 1:no_pds
%             
%             pct_bp_high_for_test(:, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*epoch_secs, 1);
%             
%         end
%         
%         if no_pds == 2 % for comp = 1:no_comparisons
%         
%             p_val(b, ch) = ranksum(pct_bp_high_for_test(:, 1), pct_bp_high_for_test(:, 1), 'tail', 'left');
%             
%         end
%         
%         subplot(1, no_chans, ch) % subplot(no_bands, 2, (b - 1)*2 + ch)
%         
%         range = [all_dimensions(@min, pct_bp_high_for_test) all_dimensions(@max, pct_bp_high_for_test)];
%         
%         range = max(eps, range);
%         
%         [bin_edges, bin_centers] = make_bins(range(1), range(2), 100, 'log'); % 25, ''); % 
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
%         title({[chan_labels{ch}, ','];[num2str(epoch_secs/60), ' Minutes of Densest High Power'];...
%             [band_labels{b}, ', p-value = ', num2str(p_val(b, ch)), ' (ranksum test)']})
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
%     save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{b}, '_hist'])
%     
% end