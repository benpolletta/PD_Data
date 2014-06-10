function PD_beta_epochs_roseplot_by_datapoint_relative_infusion

load('subjects.mat')

sampling_freq = 1000;
  
freq_labels = {'Striatal','Motor Ctx.'};

pd_labels = {'Pre-Infusion','Post-Infusion','Recovery'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    periods = sampling_freq*[0 basetime; basetime basetime + 60*20; basetime + 60*15 inf];
    
    no_periods = size(periods,1);
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    for ch = 1:2
         
        beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
        
        % Loading block numbers, epoch start and end indices.
        [blocks, epoch_starts, epoch_ends] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%*[^\n]');
        
        no_blocks = max(blocks);
        
        beta_pbf_name = [beta_listname(1:end-5),'_pbf_dp.txt'];
        
        block_starts = zeros(no_blocks,1);
        
        block_ends = zeros(no_blocks,1);
        
        datapoints = [];
        
        for b = 1:no_blocks
            
            block_indicator = blocks == b;
            
            block_start = min(epoch_starts(block_indicator));
            
            block_end = max(epoch_ends(block_indicator));
            
            block_datapoints = block_start:block_end;
            
            datapoints((end + 1):(end + length(block_datapoints))) = block_datapoints;
            
        end
        
        period_vec = zeros(size(datapoints));
        
        for i = 1:no_periods
           
            period_indicator = datapoints >= periods(i,1) & datapoints < periods(i,2); 
            
            period_vec = period_vec + period_indicator*i;
            
        end
        
        all_beta_data = load(beta_pbf_name);
        
        all_Fs = all_beta_data(:,2:3);
        
        all_Pds = all_beta_data(:,4);
        
        figure;
        
        for ch1 = 1:2
            
            for pd = 1:size(periods,1)
                
                subplot(no_periods, 2, (pd - 1)*2 + ch1)
                
                rose_plot(all_Pds(period_vec == pd), all_Fs(period_vec == pd, ch1), 20, 4:4:36);
                
                title({[folder,' ',freq_labels{ch},' High Beta Blocks, ',pd_labels{pd}];['Phase Lag by ',freq_labels{ch1},' Freq.']})
                
            end
            
        end
        
        save_as_pdf(gcf,[beta_listname(1:end-5),'_rose_dp_rel_inf'])
        
    end
          
end