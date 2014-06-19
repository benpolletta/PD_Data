function PD_beta_epochs_rel_infusion_roseplot_by_datapoint(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subject_mat)

sampling_freq = 1000;

smooth_winsize = 50;

f_bins = 8:4:32;

f_labels = textscan(num2str(f_bins), '%s', 'delimiter', ' ');
f_labels = cellstr(f_labels{1});
f_labels = f_labels(1:2:end);

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

chan_labels = {chan_labels{:}, 'Both'};

ch_label = {'ch1', 'ch2', 'ch1_ch2'};

fid_mat = nan(3,1);

all_beta_name = cell(3,1);

for ch = 1:3
    
    all_beta_name{ch} = [subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_beta_',ch_label{ch}];
    
    fid_mat(ch) = fopen([all_beta_name{ch},'_pbf_dp.txt'], 'w');
    
end

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};  
    
    base_index = basetimes(fo)*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    periods = [1 base_index; (base_index + 1) size(A,1)];
    
    no_periods = size(periods,1);
    
    for ch = 1:3
        
        beta_name = [subj_name,'_',ch_label{ch},'_beta'];
        
        figure;
        
        for pd = 1:no_periods
            
            beta_listname = [beta_name,'_',pd_label{pd},'_',par_name,'.list'];
            
            % Loading block numbers, epoch start and end indices.
            [blocks, beta_starts, beta_ends] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%*[^\n]');
            
            beta_pbf_name = [beta_listname(1:end-5),'_pbf_dp.txt'];
            
            if isempty(dir(beta_pbf_name))
                
                fid = fopen(beta_pbf_name, 'w');
                
                no_blocks = length(blocks);
                
                %% Computing smoothed frequency and phase difference for each block.
                
                for b = 1:no_blocks
                    
                    P_block = unwrap(P(beta_starts:beta_ends,:,3));
                    
                    P_diff = -diff(P_block,[],2);
                    
                    P_diff = angle(exp(sqrt(-1)*P_diff));
                    
                    Pd_flipped = [flipud(P_diff(1:smooth_winsize)); P_diff; flipud(P_diff((end-smooth_winsize+1):end))];
                    
                    Pd_conv = conv(Pd_flipped,hann(smooth_winsize)/sum(hann(smooth_winsize)),'same');
                    
                    Pd_smooth = Pd_conv((smooth_winsize + 1):(end - smooth_winsize));
                    
                    F = diff(P_block)/(2*pi*(1/sampling_freq));
                    
                    F_smooth = nan(size(Pd_smooth,1),2);
                    
                    for ch1 = 1:2
                        
                        F_flipped = [flipud(F(1:smooth_winsize,ch1)); F(:,ch1); flipud(F((end-smooth_winsize+1):end,ch1))];
                        
                        % F_conv = conv(F_flipped,hann(smooth_winsize)/sum(hann(smooth_winsize)),'same');
                        F_conv = conv(F_flipped,ones(smooth_winsize,1)/smooth_winsize,'same');
                        
                        F_smooth(:,ch1) = F_conv(smooth_winsize + (1:size(Pd_smooth,1)));
                        
                    end
                    
                    fprintf(fid, '%f\t%f\t%f\t%f\n', [b*ones(size(Pd_smooth)) F_smooth Pd_smooth]');
                    
                end
                
                fclose(fid);
                
            end
            
            all_beta_data = load(beta_pbf_name);
                
            fprintf(fid_mat(ch), '%f\t%f\t%f\t%f\t%f\n', [pd*ones(size(all_beta_data, 1), 1) all_beta_data]');
            
            if ~isempty(all_beta_data)
                
                all_Fs = all_beta_data(:,2:3);
                
                all_Pds = all_beta_data(:,4);
                
            else
                
                all_Fs = [nan nan];
                
                all_Pds = nan;
                
            end
            
            for ch1 = 1:2
                
                subplot(2, 2, (pd - 1)*2 + ch1)
                
                rose_plot(all_Pds, all_Fs(:,ch1), 20, f_bins);
                
                title({[folder,' ',chan_labels{ch},' High Beta Blocks, ',period_label{pd}];['Phase Lag by ',chan_labels{ch1},' Freq.']})
                
            end
            
        end
        
        save_as_pdf(gcf,[beta_name,'_ri_rose_dp'])
        
    end
    
end

close('all')

for ch = 1:3
    
    fclose(fid_mat(ch));
    
end

% figure(1)
% 
% index = 2;
% 
% for ch = 1:3
%             
%     all_beta_data = load([all_beta_name{ch},'_pbf_dp.txt']);
%     
%     all_pd_index = all_beta_data(:,2);
%     
%     all_Fs = all_beta_data(:,3:4);
%     
%     all_Pds = all_beta_data(:,5);
%            
%     for ch1 = 1:2
%         
%         figure(index)
%         
%         for pd = 1:length(pd_label)
%             
%             figure(index)
%             
%             subplot(1, 2, pd)
%             
%             rose_plot(all_Pds(all_pd_index == pd), all_Fs(all_pd_index == pd, ch1), 20, f_bins);
%             
%             title({[chan_labels{ch},' High Beta Blocks, ',period_label{pd}];['Phase Lag by ',chan_labels{ch1},' Freq.']})
%             
%             figure(1)
%             
%             subplot(3, 4, (ch-1)*(2 + length(pd_label)) + (pd-1)*2 + ch1)
%             
%             rose_plot(all_Pds(all_pd_index == pd), all_Fs(all_pd_index == pd, ch1), 20, f_bins);
%             
%             title({[chan_labels{ch},' High Beta Blocks, ',period_label{pd}];['Phase Lag by ',chan_labels{ch1},' Freq.']})
%             
%         end
%         
%         save_as_pdf(index, [subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_',ch_label{ch},'_by_ch',num2str(ch1),'_beta_ri_rose_dp'])
%         
%         index = index + 1;
%         
%     end
%     
% end
% 
% save_as_pdf(gcf,[subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_beta_ri_rose_dp'])
% 
% end
% 
% function F_c = categorize_freq(F, f_bins)
%     
%     [r, c] = size(F);
%     
%     F_c = zeros(r, c);
%     
%     no_f_bins = length(f_bins) - 1;
%     
%     for col = 1:c
%         
%         for f = 1:no_f_bins
%             
%             F_bin = F(:, col) >= f_bins(f) & F(:, col) < f_bins(f + 1);
%             
%             F_c(:, col) = F_c(:, col) + f*F_bin;
%             
%         end
%         
%     end
%     
% end