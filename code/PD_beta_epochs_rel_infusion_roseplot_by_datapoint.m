function PD_beta_epochs_rel_infusion_roseplot_by_datapoint

% load('initial_subjects.mat')

load('st_m1_subjects.mat')

sampling_freq = 1000;

freq_labels = {'Striatal','Motor Ctx.'};

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

for fo = 5:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);    
    
    base_index = basetimes(fo)*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    periods = [1 base_index; (base_index + 1) size(A,1)];
    
    no_periods = size(periods,1);
    
    for ch = 1:2
        
        beta_name = [subj_name,'_ch',num2str(ch),'_beta'];
        
        figure;
        
        for pd = 1:no_periods
            
            beta_listname = [beta_name,'_',pd_label{pd},'.list'];
            
            % Loading block numbers, epoch start and end indices.
            [blocks, beta_starts, beta_ends] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%*[^\n]');
            
            beta_pbf_name = [beta_listname(1:end-5),'_pbf_dp.txt'];
            
            fid = fopen(beta_pbf_name, 'w');
            
            no_blocks = length(blocks);
            
            %% Computing smoothed frequency and phase difference for each block.
            
            for b = 1:no_blocks
                
                P_block = unwrap(P(beta_starts:beta_ends,:,3));
                
                P_diff = -diff(P_block,[],2);
                
                P_diff = angle(exp(sqrt(-1)*P_diff));
                
                Pd_flipped = [flipud(P_diff(1:100)); P_diff; flipud(P_diff((end-100+1):end))];
                
                Pd_conv = conv(Pd_flipped,ones(100,1)/100,'same');
                
                Pd_smooth = Pd_conv(101:(end-100));
                
                F = diff(P_block)/(2*pi*(1/sampling_freq));
                
                F_smooth = nan(size(Pd_smooth,1),2);
                
                for ch1 = 1:2
                    
                    F_flipped = [flipud(F(1:100,ch1)); F(:,ch1); flipud(F((end-100+1):end,ch1))];
                    
                    F_conv = conv(F_flipped,ones(100,1)/100,'same');
                    
                    F_smooth(:,ch1) = F_conv(100 + (1:size(Pd_smooth,1)));
                    
                end
                
                fprintf(fid, '%f\t%f\t%f\t%f\n', [b*ones(size(Pd_smooth)) F_smooth Pd_smooth]');
                
            end
            
            fclose(fid);
            
            all_beta_data = load(beta_pbf_name);
            
            all_Fs = all_beta_data(:,2:3);
            
            all_Pds = all_beta_data(:,4);
            
            for ch1 = 1:2
                
                subplot(2, 2, (pd - 1)*2 + ch1)
                
                rose_plot(all_Pds, all_Fs(:,ch1), 20, 4:4:36);
                
                title({[folder,' ',freq_labels{ch},' High Beta Blocks, ',period_label{pd}];['Phase Lag by ',freq_labels{ch1},' Freq.']})
                
            end
            
        end
        
        save_as_pdf(gcf,[beta_name,'_ri_rose_dp'])
        
    end
    
end