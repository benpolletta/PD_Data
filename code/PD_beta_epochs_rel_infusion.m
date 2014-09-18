function PD_beta_epochs_rel_infusion(filenames, sampling_freq, chan_labels, infusion_times, outlier_lim, sd_lim, win_size, smooth_size)

% SAMPLE CALL:
% MM_beta_epochs_rel_infusion({'file1.txt','file2.txt'},1000,{'Striatum','Motor
% Ctx.'},[3000,5000],7,2,333,20000)
% 
% 'filenames' is a cell of strings, which are the filenames of files 
% containing data for picking beta segments. The data should contain two
% channels, as columns.
% 'sampling_freq' is the sampling frequency of the data.
% 'chan_labels' is a cell of strings, which are the labels of channels
% inside each file containing data.
% 'infusion_times' is a vector, with length the same as 'filenames', 
% containing the times (in datapoints) at which each data file switches
% from control to Parkinsonian behavior (or the time of carbachol
% infusion).
% 'outlier_lim' is the number of standard deviations beyond which a spike
% in the LFP is considered an outlier.
% 'sd_lim' is the number of standard deviations defining the cutoff of high
% beta power.
% 'win_size' is the minimum duration (in datapoints, so s*sampling_freq) 
% for which beta must be elevated above the cutoff, to be considered a 
% high beta segment.
% 'smooth_size' is the length of time (in datapoints, so s*sampling_freq)
% over which beta power is smoothed before applying the cutoff.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating name to call the analysis by - either name of single file in
% 'filenames', or 'All'. NB: will be over-written each time you run a
% multi-file analysis, so should be renamed after running.

if length(filenames) == 1
    
    file_label = filenames{1};
    
else

    file_label = 'All';
    
end

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

pd_label = {'pre', 'post'};

period_label = {'Pre-Infusion', 'Post-Infusion'};

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_index = {1, 2, 1:2, 1:2, 1, 2, 1, 2};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

no_b_blocks = nan(length(filenames), no_channels);

no_dps = nan(length(filenames), no_channels);

% Looping over files.

for fi = 1:length(filenames)
    
    filename = filenames{fi};
    
    infusion_index = infusion_times(fi);
    
    % Loading data and results of 'MM_bandpass'.
        
    data = load(filename);
    
    if isstruct(data)
        
        field_names = fieldnames(data);
       
        if length(field_names) ~= 1
            
            display('Input .mat file must contain one variable - data in 2 columns.')
            
            return
            
        else
            
            data = getfield(data, field_names{1});
            
        end
        
    end
        
    load([filename,'_HAP.mat'])
    
    t = (1:size(A,1))/sampling_freq;
    
    beta_amp = A(:,:,3); ba_smooth = nan(size(beta_amp)); beta_high = nan(size(beta_amp)); beta_low = nan(size(beta_amp));
    
    pd_limits = [1 infusion_index; (infusion_index + 1) min(length(t), infusion_index + 1500*sampling_freq)];
    
    % beta_cutoff = mean(beta_amp(1:base_index, :)) + sd_lim*std(beta_amp(1:base_index, :));
    
    norm_data = nan(size(data));
   
    for ch = 1:2
        
        norm_data(:,ch) = (data(:, ch) - mean(data(:, ch)))/std(data(:, ch));
        
        ba_flipped = [flipud(beta_amp(1:smooth_size, ch)); beta_amp(:,ch); flipud(beta_amp((end-smooth_size+1):end,ch))];
        
        ba_conv = conv(ba_flipped, hann(smooth_size)/sum(hann(smooth_size)), 'same');
        
        ba_smooth(:, ch) = ba_conv((smooth_size+1):(end-smooth_size));
        
        beta_h_cutoff = mean(ba_smooth(1:infusion_index, ch)) + sd_lim*std(ba_smooth(1:infusion_index, ch));
        
        beta_l_cutoff = mean(ba_smooth(1:infusion_index, ch)) - sd_lim*std(ba_smooth(1:infusion_index, ch));
        
        % beta_h_cutoff = quantile(ba_smooth(1:base_index, ch), sd_lim);
        % 
        % beta_l_cutoff = quantile(ba_smooth(1:base_index, ch), 1 - sd_lim);
        
        beta_high(:, ch) = ba_smooth(:, ch) >= beta_h_cutoff & abs(norm_data(:, ch)) < outlier_lim;
        
        beta_low(:, ch) = ba_smooth(:, ch) <= beta_l_cutoff & abs(norm_data(:, ch)) < outlier_lim;
        
    end
    
    beta_high(:, 3) = sum(beta_high, 2) >= 2;
    
    beta_high(:, 4) = sum(beta_high, 2) >= 1;
    
    beta_high(:, 5) = beta_high(:, 1) & ~beta_high(:, 2);
    
    beta_high(:, 6) = ~beta_high(:, 1) & beta_high(:, 2);
    
    beta_high(:, 7) = beta_high(:, 1) & beta_low(:, 2);
    
    beta_high(:, 8) = beta_low(:, 1) & beta_high(:, 2);
    
    figure;
    
    for ch = 1:no_channels
        
        for pd = 1:size(pd_limits,1)
            
            beta_listname = [filename,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,'.list'];
            
            fid_list = fopen(beta_listname, 'w');
            
            fid_A_list = fopen([beta_listname(1:end-5), '_A.list'],'w');
            
            fid_P_list = fopen([beta_listname(1:end-5), '_P.list'],'w');
            
            fid_win_list = fopen([beta_listname(1:end-5), '_win.list'],'w');
            
            subplot(no_channels + 1, 2, (ch-1)*2 + pd)
            
            plot(t, beta_amp(:, ch_index{ch}), 'k', t, ba_smooth(:, ch_index{ch}), 'b')
            
            % plot(t(pd_limits(pd,1):pd_limits(pd,2)), beta_amp(pd_limits(pd,1):pd_limits(pd,2), ch_index{ch}), 'k',...
            %     t(pd_limits(pd,1):pd_limits(pd,2)), ba_smooth(pd_limits(pd,1):pd_limits(pd,2), ch_index{ch}), 'b')
            
            title([filename, ' ', chan_labels{ch}, ' Beta Segments ', period_label{pd}])
            
            axis tight
            
            hold on
            
            % ba_pd = ba_smooth(pd_limits(pd,1):pd_limits(pd,2));
            
            % beta_cutoff = mean(ba_pd) + sd_lim*std(ba_pd);
            
            % beta_high = ba_pd >= beta_cutoff;
            
            bh_pd = beta_high(pd_limits(pd,1):pd_limits(pd,2), ch);
            
            if any(bh_pd == 0)
                
                dbh = diff(bh_pd);%beta_high);
                
                beta_start = find(dbh == 1) + 1 + pd_limits(pd,1) - 1;
                
                beta_end = find(dbh == -1) + pd_limits(pd,1) - 1;
                
            else
                
                beta_start = pd_limits(pd, 1);
                
                beta_end = pd_limits(pd, 2);
                
            end
            
            if ~isempty(beta_end) || ~isempty(beta_start)
                
                if isempty(beta_start)
                    
                    beta_start = [pd_limits(pd,1); beta_start];
                    
                end
                
                if isempty(beta_end)
                    
                    beta_end = [beta_end; pd_limits(pd, 2)];
                    
                end
                
                if beta_end(1) < beta_start(1)
                    
                    beta_start = [pd_limits(pd,1); beta_start];
                        
                end
                
                if beta_start(end) > beta_end(end)
                    
                    beta_end = [beta_end; pd_limits(pd,2)];
                    
                end
                
                beta_end(end) = min(beta_end(end), length(t));
                
                beta_blocks = [beta_start beta_end];
                
                beta_lengths = diff(beta_blocks, [], 2) + 1;
                
                beta_blocks(beta_lengths < win_size, :) = [];
                
                beta_lengths = diff(beta_blocks, [], 2) + 1;
                
                for b = 1:size(beta_blocks,1)
                        
                    beta_name = [filename,'_',par_name,'_ch',num2str(ch),'_beta_',pd_label{pd},'_block',num2str(b),'.txt'];
                        
                    plot(t(beta_blocks(b,1):beta_blocks(b,2)), beta_amp(beta_blocks(b,1):beta_blocks(b,2), ch_index{ch}), 'g')

                    plot(t(beta_blocks(b,1):beta_blocks(b,2)), ba_smooth(beta_blocks(b,1):beta_blocks(b,2), ch_index{ch}), 'r')
                    
                    no_epochs = floor(beta_lengths(b)/win_size);
                    
                    for e = 1:no_epochs
                        
                        epoch_name = [beta_name,'_epoch',num2str(e),'.txt'];
                        
                        A_name = [beta_name,'_epoch',num2str(e),'_A.txt'];
                        
                        P_name = [beta_name,'_epoch',num2str(e),'_P.txt'];
                        
                        epoch_start = beta_blocks(b,1) + (e-1)*win_size;
                        
                        epoch_end = beta_blocks(b,1) + e*win_size - 1;
                        
                        fid = fopen(epoch_name, 'w');
                        
                        fprintf(fid, '%f\t%f\n', data(epoch_start:epoch_end, :)');
                        
                        fclose(fid);
                        
                        fid = fopen(A_name, 'w');
                        
                        fprintf(fid, '%f\t%f\n', A(epoch_start:epoch_end, :, 3)');
                        
                        fclose(fid);
                        
                        fid = fopen(P_name, 'w');
                        
                        fprintf(fid, '%f\t%f\n', P(epoch_start:epoch_end, :, 3)');
                        
                        fclose(fid);
                        
                        fprintf(fid_list, '%s\n', epoch_name);
                        
                        fprintf(fid_A_list, '%s\n', A_name);
                        
                        fprintf(fid_P_list, '%s\n', P_name);
                        
                        fprintf(fid_win_list, '%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n', b, beta_blocks(b,1), beta_blocks(b,2), e, epoch_start, epoch_end, median(beta_amp(epoch_start:epoch_end,:)));
                        
                    end
                    
                end
                
                no_b_blocks(fi, ch, pd) = size(beta_blocks, 1);
                
                beta_lengths = diff(beta_blocks, [], 2) + 1;
                
                no_dps(fi, ch, pd) = sum(beta_lengths);
                
            end
            
        end
        
    end
    
    subplot(no_channels + 1, 2, 15)
    
    bar(reshape(no_b_blocks(fi, :, :), no_channels, 2))
    
    title([filename, ' Number Beta Segments'])
    
    set(gca,'XTickLabel',ch_label)
    
    legend(period_label)
    
    subplot(no_channels + 1, 2, 16)
    
    bar(reshape(no_dps(fi, :)/sampling_freq, no_channels, 2))
    
    title([filename, ' Length Beta Segments (s)'])
    
    set(gca,'XTickLabel',ch_label)
    
    legend(period_label)
    
    save_as_pdf(gcf,[filename, '_beta_', par_name])
          
end

figure

subplot(1, 2, 1)

if size(no_b_blocks, 1) > 1
    
    no_b_blocks = sum(no_b_blocks);
    
end

bar(reshape(no_b_blocks, no_channels, 2))
    
title('Total Number Beta Segments')

set(gca,'XTickLabel',ch_label)

legend(period_label)

subplot(1, 2, 2)

if size(no_dps, 1) > 1
    
    no_dps = sum(no_b_blocks);
    
end

bar(reshape(no_dps/sampling_freq, no_channels, 2))

title('Total Length Beta Segments (s)')

set(gca,'XTickLabel',ch_label)

legend(period_label)

save_as_pdf(gcf, [file_label, '_beta_', par_name])
