name = 'STR_w_M1_1-200Hz_3-21cycles_7bands_kmeans_pct_power_post_vs_pre_summary';

st_m1 = power_post_vs_pre_summary('st_m1_subjects.mat', '_kmeans', 150, '_by_STR', '_pct', freqs, no_cycles, bands, 4);
st_m1_ali = power_post_vs_pre_summary('st_m1_ali_subjects.mat', '_kmeans', 150, '_by_STR', '_pct', freqs, no_cycles, bands, 4);
st_m1_ali2 = power_post_vs_pre_summary('st_m1_ali2_subjects.mat', '_kmeans', 150, '_by_STR', '_pct', freqs, no_cycles, bands, 4);

for stat = 1:3, st_m1_ali(1:2, :, stat) = flipud(st_m1_ali(1:2, :, stat)); end
for stat = 1:3, st_m1_ali2(1:2, :, stat) = flipud(st_m1_ali2(1:2, :, stat)); end

increase_summary = st_m1;
for stat = 1:3, increase_summary(:, 10:12, stat) = st_m1_ali(:, :, stat); end
for stat = 1:3, increase_summary(:, 13:14, stat) = st_m1_ali2(:, :, stat); end

save([name, '.mat'], 'increase_summary')

stat_labels = {'_total', '_start', '_end'}; no_stats = length(stat_labels);

for stat = 1:no_stats
   
    fid = fopen([name, stat_labels{stat}, '.txt'], 'w');
    
    fprintf(fid, make_format(3, 'f'), increase_summary(:, :, stat));
    
    fclose(fid);
    
end

no_blocks = 3;

stat_labels = {'Total Duration', 'Start Time', 'End Time'};

chan_labels = {'Striatum', 'M1'};

figure

for stat = 1:no_stats
    
    subplot(2, no_stats, stat)

    plot((1:no_blocks)', increase_summary(:, :, stat))
    
    hold on
    
    for block = 1:no_blocks
    
        plot(block, increase_summary(block, :, stat)', 'o')
        
    end
    
    set(gca, 'XLim', [.5 no_stats + .5], 'XTick', 1:no_stats, 'XTickLabel', {chan_labels{:}, 'Overlap'}, 'FontSize', 14)
    
    ylabel('Time (min.)')
    
    title(stat_labels{stat})
    
    subplot(2, no_stats, no_stats + stat)
    
    barwitherr(nanstd(increase_summary(:, :, stat)'), nanmean(increase_summary(:, :, stat)'));
    
    set(gca, 'XTickLabel', {chan_labels{:}, 'Overlap'}, 'FontSize', 14)
    
    ylabel('Time (min.)')
    
end

save_as_pdf(gcf, 'STR_w_M1_1-200Hz_3-21cycles_7bands_kmeans_pct_power_post_vs_pre_summary')