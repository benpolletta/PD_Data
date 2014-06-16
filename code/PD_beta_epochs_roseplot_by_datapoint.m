function PD_beta_epochs_roseplot_by_datapoint(subject_mat)

load(subject_mat)

sampling_freq = 1000;

f_bins = 8:4:32;

fid_vec = nan(2,1); all_beta_name = cell(2,1);

for ch = 1:2
    
    all_beta_name{ch} = [subject_mat(1:(end-length('_subjects.mat'))),'_beta_ch',num2str(ch)];
    
    fid_vec(ch) = fopen([all_beta_name{ch},'_pbf_dp.txt'], 'w');
    
end

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    for ch = 1:2
         
        beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
        
        % Loading block numbers, epoch start and end indices.
        [blocks, epoch_starts, epoch_ends] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%*[^\n]');
        
        beta_pbf_name = [beta_listname(1:end-5),'_pbf_dp.txt'];
        
        if isempty(dir(beta_pbf_name))
            
            fid = fopen(beta_pbf_name,'w');
            
            no_blocks = max(blocks);
            
            for b = 1:no_blocks
                
                figure;
                
                block_indicator = blocks == b;
                
                block_start = min(epoch_starts(block_indicator));
                
                block_end = max(epoch_ends(block_indicator));
                
                %% Plotting phase diff. and frequency.
                
                t = (block_start:block_end)'/sampling_freq;
                
                P_block = unwrap(P(block_start:block_end,:,3));
                
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
                
                subplot(3,1,1)
                
                [ax, ~, h2] = plotyy(t,F_smooth,t,Pd_smooth/pi);
                
                title([folder,' Channel ',num2str(ch),' Block ',num2str(b),', Beta Frequency and Phase Lag'])
                
                xlabel('Time (s)')
                
                axis(ax(1),'tight'), axis(ax(2),'tight')
                
                set(ax,{'ycolor'},{'black';'black'})
                
                set(get(ax(1),'YLabel'),'String','\beta Frequency (Hz)')
                
                legend(ax(1), {'Striatum','Motor Ctx.'})
                
                set(get(ax(2),'YLabel'),'String','Phase Diff. (Motor - Striatal, \pi)')
                
                set(h2,'Color','k')
                
                for ch1 = 1:2
                    
                    subplot(2, 2, 2 + ch1)
                    
                    rose_plot(P_diff, F_smooth(:,ch1), 20, f_bins);
                    
                    title({[folder,' Block ',num2str(b)];['Phase Lag by ',chan_labels{ch1},' Freq.']})
                    
                end
                
                save_as_pdf(gcf,[beta_listname(1:end-5),'_block',num2str(b),'_rose_dp'])
                
                close(gcf)
                
            end
            
            fclose(fid);
            
        end
        
        all_beta_data = load(beta_pbf_name);
        
        fprintf(fid_vec(ch), '%f\t%f\t%f\t%f\n', all_beta_data');
        
        all_Fs = all_beta_data(:,2:3);
        
        all_Pds = all_beta_data(:,4);
        
        figure;
        
        for ch1 = 1:2
        
            subplot(1, 2, ch1)
            
            rose_plot(all_Pds, all_Fs(:,ch1), 20, f_bins);
        
            title({[folder,' ',chan_labels{ch},' High Beta Blocks'];['Phase Lag by ',chan_labels{ch1},' Freq.']})
           
        end
        
        save_as_pdf(gcf,[beta_listname(1:end-5),'_rose_dp'])
        
    end
          
end

for ch = 1:2
    
    fclose(fid_vec(ch));
    
end

figure;

for ch = 1:2
    
    all_beta_data = load([all_beta_name{ch},'_pbf_dp.txt']);
    
    all_Fs = all_beta_data(:,2:3);
    
    all_Pds = all_beta_data(:,4);
    
    for ch1 = 1:2
        
        subplot(2, 2, (ch-1)*2 + ch1)
        
        rose_plot(all_Pds, all_Fs(:,ch1), 20, f_bins);
        
        title({[chan_labels{ch},' High Beta Blocks'];['Phase Lag by ',chan_labels{ch1},' Freq.']})
        
    end
    
end

save_as_pdf(gcf,[subject_mat(1:(end-length('_subjects.mat'))),'_beta_rose_dp'])