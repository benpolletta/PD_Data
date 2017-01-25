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

stat_labels = {'_total', '_start', '_end'};

for stat = 1:3
   
    fid = fopen([name, stat_labels{stat}, '.txt'], 'w');
    
    fprintf(fid, make_format(3, 'f'), increase_summary(:, :, stat));
    
    fclose(fid);
    
end