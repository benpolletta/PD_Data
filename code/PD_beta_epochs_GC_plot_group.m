function PD_beta_epochs_GC_plot_group(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

period_label = {'Pre-Injection','Post-Injection'};
pd_label = {'pre', 'post'};

list_suffix = {'', '_A'};

no_lists = length(list_suffix);

f = 1000*(0:win_size)/win_size;

f_indices{1} = f <= 31 & f >= 9; %f <= 100
f_indices{2} = f <= 100; %40;

f_length = win_size + 2;

measure_label = {'',', p-Value'};
no_measures = length(measure_label);

list_label = {'GC of LFP','GC of Beta Amp.'};

pd_colors = {'b','r','g','k'};

dir_labels = {'ch1toch2', 'ch2toch1'};
direction_labels{1} = [chan_labels{1}, ' --> ', chan_labels{2}];
direction_labels{2} = [chan_labels{2}, ' --> ', chan_labels{1}];
no_directions = length(direction_labels);

pd_dir_label = cell(4, 1);
for pd = 1:2
    for dir = 1:2
        pd_dir_label{2*pd - (2 - dir)} = [direction_labels{dir}, ', ', period_label{pd}];
    end
end
        
All_mean_GC_spec = nan(f_length, 2, no_channels, no_lists);

All_std_GC_spec = nan(f_length, 2, 2, no_channels, no_lists);

for ch = 1:4%no_channels
    
    for list = 1:no_lists
        
        figure
        
        mean_GC_spec = nan(f_length, 2);
        
        std_GC_spec = nan(f_length, 1, 2);
        
        for pd = 1:2
            
            listname = ['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_', ch_label{ch}, '_', pd_label{pd}, list_suffix{list}, '_GC'];
            
            All_GC = load([listname, '.txt']);
            
            if isempty(All_GC)
            
                All_GC = nan(1, 6);
                
            end
                
            All_GC_spec = nan(size(All_GC, 1), f_length, 2);
            
            for dir = 1:no_directions
                
                dir_GC_spec = load([listname, 'spec_', dir_labels{dir}, '.txt']);
                
                if ~isempty(dir_GC_spec)
                    
                    All_GC_spec(:, :, dir) = dir_GC_spec;
                    
                end
                
            end
            
            errors = All_GC(:, 2);
            
            All_GC(errors == 1, :) = [];
            
            All_GC_spec(errors == 1, :, :) = [];
            
            no_epochs = size(All_GC, 1);
            
            for measure = 1:no_measures
                
                for dir = 1:no_directions
                    
                    subplot(2, 2, measure)
                    
                    [histogram, bins] = hist(All_GC(:, 2 + 2*measure - (2 - dir)), 50);
                    
                    plot(bins, histogram/sum(histogram), pd_colors{2*pd - (2 - dir)})
                    
                    hold on
                    
                    title({[chan_labels{ch}, ' High Beta Blocks'];['Histogram, ', list_label{list}, measure_label{measure}]})
                    
                    if pd == 2 && dir == 2 && measure == 1
                        
                        legend(pd_dir_label)
                        
                    end
                    
                end
                
            end
                
            All_GC_spec_diff = -diff(All_GC_spec, [], 3);
            
            if size(All_GC_spec_diff, 1) > 1
                
                mean_GC_spec(:, pd) = reshape(nanmean(All_GC_spec_diff), f_length, 1);
                
                std_GC_spec(:, 1, pd) = reshape(nanstd(All_GC_spec_diff)/sqrt(no_epochs), f_length, 1);
                
            else
               
                mean_GC_spec(:, pd) = reshape(All_GC_spec_diff, f_length, 1);
                
                std_GC_spec(:, 1, pd) = nan(f_length, 1);
                
            end
            
        end
        
        All_mean_GC_spec(:, :, ch, list) = mean_GC_spec;
        
        std_GC_spec = repmat(std_GC_spec, [1 2 1]);
        
        All_std_GC_spec(:, :, :, ch, list) = std_GC_spec;
        
        subplot(2, 1, 2)
        
        boundedline(f(f_indices{list}), mean_GC_spec(f_indices{list}, :), std_GC_spec(f_indices{list}, :, :))
        
        legend(period_label)
        
        axis tight
        
        title({[chan_labels{ch}, ' High Beta Blocks'];['Spectral ', list_label{list}, ', (', direction_labels{1}, ') - (', direction_labels{2}, ')']})
        
        save_as_pdf(gcf, listname)
        
    end
    
    close('all')
    
end

save(['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_', par_name, '_GC_spec_stats.mat'], 'All_mean_GC_spec', 'All_std_GC_spec')