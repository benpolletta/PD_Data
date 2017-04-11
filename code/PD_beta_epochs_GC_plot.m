function PD_beta_epochs_GC_plot(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_index = {1, 2, 1:2, 1:2, 1, 2, 1, 2};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

period_label = {'Pre-Injection','Post-Injection'};
pd_label = {'pre', 'post'};

list_suffix = {'', '_A'};

no_lists = length(list_suffix);

f = 1000*(0:win_size)/win_size;

f_indices{1} = f <= 35 & f >= 5;
f_indices{2} = f <= 40;

f_length = win_size + 2;

measure_label = {'',', p-Value'};
no_measures = length(measure_label);

list_label = {'GC of LFP','GC of Beta Amp.'};

pd_colors = {'b','r','g','k'};

direction_labels{1} = [chan_labels{1}, ' --> ', chan_labels{2}];
direction_labels{2} = [chan_labels{2}, ' --> ', chan_labels{1}];
no_directions = length(direction_labels);

for pd = 1:2
    for dir = 1:2
        pd_dir_label{2*pd - (2 - dir)} = [direction_labels{dir}, ', ', period_label{pd}];
    end
end

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
        
    for ch = 1:no_channels
        
        for list = 1:no_lists
            
            figure;
            
            for pd = 1:2
                
                beta_listname = [subj_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,list_suffix{list},'.list'];
                % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
                
                load([beta_listname(1:end-5),'_GC.mat'],'All_GC','All_GC_spec')
                
                no_epochs = size(All_GC, 1);
                
                errors = All_GC(:, 2);
                
                All_GC(errors == 1, :) = [];
                
                for measure = 1:no_measures
                    
                    for dir = 1:no_directions
                
                        subplot(2, 2, measure)
                        
                        [histogram, bins] = hist(All_GC(:, 2 + 2*measure - (2 - dir)), 50);
                        
                        plot(bins, histogram/sum(histogram), pd_colors{2*pd - (2 - dir)})
                        
                        hold on
                        
                        title({[folder, ', ', chan_labels{ch}, ' High Beta Blocks'];['Histogram, ', list_label{list}, measure_label{measure}]})
                        
                        if pd == 2 && dir == 2 && measure == 1
                           
                            legend(pd_dir_label)
                            
                        end
                        
                    end
                    
                end
                
                mean_GC_spec = reshape(nanmean(All_GC_spec), f_length, 2);
                
                std_GC_spec = reshape(nanstd(All_GC_spec)/sqrt(no_epochs), f_length, 1, 2);
                std_GC_spec = repmat(std_GC_spec, [1 2 1]);
                
                subplot(2, 2, 2 + pd)
                
                boundedline(f(f_indices{list}), mean_GC_spec(f_indices{list}, :), std_GC_spec(f_indices{list}, :, :))
                
                legend(direction_labels)
                
                axis tight
                
                title({[folder, ', ', chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];['Spectral ', list_label{list}]})
                
            end
            
            save_as_pdf(gcf, beta_listname(1:end-5))
            
        end
        
    end
    
end