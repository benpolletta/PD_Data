function PD_beta_blocks_rel_infusion_pre_post_plot(subject_mat, epoch_secs, pd_handle, freqs, no_cycles, bands)

subj_mat_name = subject_mat(1:(end - length('_subjects.mat')));

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

close('all')

load(subject_mat)

no_folders = length(folders);

if exist([folders{1}, '/', prefixes{1}, '_wt.mat'], 'file')
    
    load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')
    
else
    
    sampling_freq = 500;
    
end

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

if ~isempty(epoch_secs)
    
    epoch_secs_label = ['_', num2str(epoch_secs/60), '_min_secs'];
    
    epochs_label = ['_', num2str(epoch_secs/60)];
    
else
    
    epoch_secs = 5*str2num(pd_handle(2:(end - length('trials'))));
    
    epoch_secs_label = '_laser';
    
    epochs_label = '';
    
end
    
load([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high', epoch_secs_label, pd_handle, '.mat'])

power_flag = exist('BP_sec', 'var') | strcmp(pd_handle((end - 5):end), '_power') | strcmp(pd_handle, '_pct_power');

if exist('BP_sec', 'var')
    
    pct_bp_high = BP_sec;
    
end

%% Group bar plots.

if no_pds > 1
    
    no_comparisons = nchoosek(no_pds, 2);
    
    comparisons = nchoosek(1:no_pds, 2);

else
    
    no_comparisons = 0;
    
end

if no_pds > 1
    
    no_comparisons = nchoosek(no_pds, 2);
    
    p_val_labels = comp_labels(pd_labels);
    
else
    
    no_comparisons = 0;
    
    p_val_labels = {};
    
end
    
p_vals = nan(no_comparisons, no_bands, no_chans, 2);

format = make_format(5*no_pds, 'f');

format = [format(1:(end - 2)), '\t', make_format(2*no_comparisons, '.3E')];

format = ['%s\t', format];

stat_labels = {'Band'};

stats = {'Mean', 'S.E.', 'Median', 'Q1', 'Q3'};

for s = 1:length(stats)
    
    for p = 1:no_pds
        
        stat_labels = {stat_labels{:}, stats{s}};
        
    end
    
end
    
fid = nan(no_chans, 1);

for ch = 1:no_chans
    
    fid(ch) = fopen([subj_mat_name, BP_suffix, epochs_label, '_mins_ranksum',...
        pd_handle, '_', chan_labels{ch}, '_stats.txt'], 'w');
    
    fprintf(fid(ch), make_format(1 + 5*no_pds + 2*no_comparisons, 's'), stat_labels{:}, p_val_labels{:});
    
    fprintf(fid(ch), make_format(1 + 5*no_pds, 's'), '', pd_labels{:}, pd_labels{:},...
        pd_labels{:}, pd_labels{:}, pd_labels{:});
    
end

for b = 1:no_bands
    
    figure
    
    for ch = 1:no_chans
        
        pct_bp_high_for_test = nan(no_folders*epoch_secs, no_pds);
        
        for pd = 1:no_pds
            
            pct_bp_high_for_test(:, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*epoch_secs, 1);
            
        end
        
        if no_pds == 2 %for comp = 1:no_comparisons
        
            p_vals(:, b, ch) = ranksum(pct_bp_high_for_test(:, 1), pct_bp_high_for_test(:, 2), 'tail', 'left');
                
        elseif no_pds > 1
            
            for comp = 1:no_comparisons
                
                if any(~isnan(pct_bp_high_for_test(:, comparisons(comp, 1)))) && any(~isnan(pct_bp_high_for_test(:, comparisons(comp, 2))))
                    
                    p_vals(comp, b, ch, 1) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'left');
                    
                    p_vals(comp, b, ch, 2) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'right');
                    
                end
                
            end
            
        end
            
        fprintf(fid(ch), format, short_band_labels{b}, nanmean(pct_bp_high_for_test), nanstd(pct_bp_high_for_test)/sqrt(epoch_secs),...
            nanmedian(pct_bp_high_for_test), quantile(pct_bp_high_for_test, .25), quantile(pct_bp_high_for_test, .75),...
            reshape(p_vals(:, b, ch, :), 1, 2*no_comparisons));
        
        subplot(1, no_chans, ch), % subplot(no_bands, 2, (b - 1)*2 + ch)
        
        if power_flag
            
            % h = boxplot(pct_bp_high_for_test);
            
            md = nanmedian(pct_bp_high_for_test);
            
            for p = 1:no_pds
                
                q1(p) = nanmedian(pct_bp_high_for_test(pct_bp_high_for_test(:, p) <= md(p), p));
                
                q3(p) = nanmedian(pct_bp_high_for_test(pct_bp_high_for_test(:, p) >= md(p), p));
                
            end
            
            if length(md) == 1
                
                md = [md nan]; q1 = [q1 nan]; q3 = [q3 nan];
                
            end
            
            err(:, :, 1) = q1 - md; err(:, :, 2) = q3 - md;
            
            h = barwitherr(err, md, 0.6);
            
            % mn = nanmean(pct_bp_high_for_test);
            % 
            % sd = nanstd(pct_bp_high_for_test)*diag(1./sqrt(sum(~isnan(pct_bp_high_for_test))));
            % 
            % if length(mn) == 1 % Stupid Matlab can't tell the difference between x & y and y & width.
            % 
            %     mn = [mn nan]; sd = [sd nan];
            % 
            % end
            % 
            % h = barwitherr(sd, mn, 0.6);
            
            ylabel('Beta Power')
            
        else
        
            pct_bp_high_for_test = max(pct_bp_high_for_test, eps);
            
            mn = nanmean(log(pct_bp_high_for_test));
            
            sd = nanstd(log(pct_bp_high_for_test))*diag(1./sqrt(sum(~isnan(log(pct_bp_high_for_test)))));
            
            if length(mn) == 1 % Stupid Matlab can't tell the difference between x & y and y & width.
                
                mn = [mn mn - 5]; sd = [sd nan];
                
            end
           
            h = barwitherr(sd, mn, 0.6, 'BaseValue', min(nanmean(log(pct_bp_high_for_test))) - 5);
        
            ylabel('Beta Density')
            
        end
        
        box off
        
        set(gca, 'XTick', 1:no_pds, 'XTickLabel', pd_labels)
        
        bar_pos = get_bar_pos(h);
        
        for comp_type = 1:2
            
            bar_pairs = {};
            
            for comp = 1:no_comparisons
                
                if p_vals(comp, b, ch, comp_type) < .05
                    
                    bar_pairs = {bar_pairs{:}, [bar_pos(comparisons(comp, 1)), bar_pos(comparisons(comp, 2))]};
                    
                end
                
            end
            
            sigstar(bar_pairs, p_vals(p_vals(:, b, ch, comp_type) < .05, b, ch, comp_type))
            
        end
        
        if no_pds == 2
            
            title({[chan_labels{ch}, ','];[num2str(epoch_secs/60), ' Minutes of Densest High Power'];...
                [band_labels{b}, ', p-value = ', num2str(p_vals(:, b, ch)), ' (ranksum test)']})
            
        else
            
            title({[chan_labels{ch}, ', ' band_labels{b}];[num2str(epoch_secs/60), ' Minutes of Densest High Power']})
            
        end
            
        xlim([0 (no_pds + 1)])
        
        set(gca, 'XTick', 1:no_pds, 'XTickLabel', pd_labels)
        
    end

    save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high', ...
        epoch_secs_label, '_mins_', short_band_labels{b}, '_barplot', pd_handle])
    
end

fclose('all');
        
save([subj_mat_name, BP_suffix, epochs_label, '_ranksum', pd_handle, '_pvals.mat'], 'p_vals')

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
%     save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{b}, '_hist', pd_handle])
%     
% end

end

function pos_bars = get_bar_pos(handle)

    for i = 1:length(handle)

        x = get(get(handle(i), 'children'), 'xdata');

        x = mean(x([1 3],:));

        pos_bars(i,:) = x;

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