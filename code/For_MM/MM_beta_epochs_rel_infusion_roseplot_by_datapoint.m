function MM_beta_epochs_rel_infusion_roseplot_by_datapoint(filenames, sampling_freq, chan_labels, injection_times, outlier_lim, sd_lim, win_size, smooth_size)

% Computes instantaneous frequency and instantaneous phase difference for
% high beta epochs.
%
% SAMPLE CALL:
% MM_beta_epochs_rel_infusion_roseplot_by_datapoint({'file1.txt','file2.txt'},1000,{'Striatum','Motor
% Ctx.'},7,2,333,20000)
% 
% INPUTS:
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

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

smooth_winsize = 50;

f_bins = 9.5:1:30.5;

% f_labels = textscan(num2str(f_bins), '%s', 'delimiter', ' ');
% f_labels = cellstr(f_labels{1});
% f_labels = f_labels(1:2:end);

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

% ch_index = {1, 2, 1:2, 1, 2, 1, 2};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

% Creating name to call the analysis by - either name of single file in
% 'filenames', or 'All'. NB: will be over-written each time you run a
% multi-file analysis, so should be renamed after running.

if length(filenames) == 1
    
    file_label = filenames{1};
    
else

    file_label = 'All';
    
end

% Name for each channel.

fid_mat = nan(no_channels,1);

all_beta_name = cell(no_channels,1);

for ch = 1:no_channels
    
    all_beta_name{ch} = [file_label,'_',par_name,'_beta_',ch_label{ch}];
    
    fid_mat(ch) = fopen([all_beta_name{ch},'_pbf_dp.txt'], 'w');
    
end

% Loop over filenames.

for fi = 1:length(filenames)
    
    filename = filenames{fi};
    
    injection_index = injection_times(fi);
    
    % Loading bandpassed data.
    
    load([filename,'_HAP.mat'])
    
    % Matrix of periods - pre-'infusion' and post-'infusion'.
    
    periods = [1 injection_index; (injection_index + 1) size(A,1)];
    
    no_periods = size(periods,1);
    
    % Loop over channels and other combinations.
    
    for ch = 1:4 %no_channels
        
        beta_name = [filename,'_',ch_label{ch},'_beta'];
        
        % figure;
        
        % Loop over pre- and post-infusion periods.
        
        for pd = 1:no_periods
            
            beta_listname = [beta_name,'_',pd_label{pd},'_',par_name,'.list'];
            
            % Loading block numbers, epoch start and end indices.
            [blocks, beta_starts, beta_ends, epochs, epoch_starts, epoch_ends] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%f%f%f%*[^\n]');
            
            beta_pbf_name = [beta_listname(1:end-5),'_pbf_dp.txt'];
            
            % if isempty(dir(beta_pbf_name))
                
                fid = fopen(beta_pbf_name, 'w');
                
                no_blocks = max(blocks);
                
                %% Computing smoothed frequency and phase difference for each block.
                
                for b = 1:no_blocks
                    
                    block_start = min(beta_starts(blocks == b));
                    
                    block_end = max(beta_ends(blocks == b));
                    
                    % Calculating smoothed phase difference between
                    % channels.
                    
                    P_block = unwrap(P(block_start:block_end, :, 3)); % Selecting block of phases.
                    
                    P_diff = -diff(P_block,[],2); % Computing phase difference.
                    
                    P_diff = angle(exp(sqrt(-1)*P_diff)); % Finding representative for phase difference in [-pi, pi].
                    
                    % Smoothing phase difference.
                    Pd_flipped = [flipud(P_diff(1:smooth_winsize)); P_diff; flipud(P_diff((end-smooth_winsize+1):end))];
                    
                    Pd_conv = angle(conv(exp(sqrt(-1)*Pd_flipped), hann(smooth_winsize)/sum(hann(smooth_winsize)), 'same'));
                    
                    Pd_smooth = Pd_conv((smooth_winsize + 1):(end - smooth_winsize));
                    
                    % Calculating smoothed instantaneous frequency of each
                    % channel.
                    
                    F = diff(P_block)/(2*pi*(1/sampling_freq)); % Calculating instantaneous frequency as difference (in time) of phases.
                    
                    % Smoothing instantaneous frequency.
                    F_smooth = nan(size(Pd_smooth,1),2);
                    
                    for ch1 = 1:2
                        
                        F_flipped = [flipud(F(1:smooth_winsize,ch1)); F(:,ch1); flipud(F((end-smooth_winsize+1):end,ch1))];
                        
                        % F_conv = conv(F_flipped,hann(smooth_winsize)/sum(hann(smooth_winsize)),'same');
                        F_conv = conv(F_flipped,ones(smooth_winsize,1)/smooth_winsize,'same');
                        
                        F_smooth(:,ch1) = F_conv(smooth_winsize + (1:size(Pd_smooth,1)));
                        
                    end
                    
                    % Printing frequency by block.
                    
                    fprintf(fid, '%f\t%f\t%f\t%f\n', [b*ones(size(Pd_smooth)) F_smooth Pd_smooth]'); % Printing smoothed phase and freq.
                    
                    % Printing frequency by epoch.
                    
                    block_name = [filename,'_',par_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_block',num2str(b)];
                    
                    no_epochs = max(epochs(blocks == b));
                    
                    block_e_starts = epoch_starts(blocks == b) - min(epoch_starts(blocks == b)) + 1;
                    
                    block_e_ends = epoch_ends(blocks == b) - min(epoch_starts(blocks == b)) + 1;
                    
                    for e = 1:no_epochs
                        
                        % display(sprintf('Block %d Epoch %d', b, e))
                        
                        e_start = block_e_starts(e);
                        
                        e_end = block_e_ends(e);
                       
                        epoch_name = [block_name,'_epoch',num2str(e),'_F.txt'];
                        
                        fid_epoch = fopen(epoch_name, 'w');
                        
                        fprintf(fid_epoch, '%f\t%f\n', F_smooth(e_start:e_end, :)');
                        
                        fclose(fid_epoch);
                        
                    end
                    
                end
                
                fclose(fid);
                
            % end
            
            %% Collecting and saving data.
            
            all_beta_data = load(beta_pbf_name);
            
            format = make_format(size(all_beta_data, 2) + 1, 'f');
                
            fprintf(fid_mat(ch), format, [pd*ones(size(all_beta_data, 1), 1) all_beta_data]');
            
            % if ~isempty(all_beta_data)
            % 
            %     all_Fs = all_beta_data(:,2:3);
            % 
            %     all_Pds = all_beta_data(:,4);
            % 
            % else
            % 
            %     all_Fs = [nan nan];
            % 
            %     all_Pds = nan;
            % 
            % end
            % 
            % for ch1 = 1:2
            % 
            %     subplot(2, 2, (pd - 1)*2 + ch1)
            % 
            %     rose_plot(all_Pds, all_Fs(:,ch1), 20, f_bins);
            % 
            %     title({[folder,' ',chan_labels{ch},' High Beta Blocks, ',period_label{pd}];['Phase Lag by ',chan_labels{ch1},' Freq.']})
            % 
            % end
            
        end
        
        % save_as_pdf(gcf,[beta_name,'_ri_rose_dp'])
        
    end
    
end

close('all')

% Closing group files collecting smoothed phase and freq.

for ch = 1:no_channels
    
    fclose(fid_mat(ch));
    
end