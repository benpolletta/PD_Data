function MM_beta_epochs_coh_mtm_collect(filenames, chan_labels, infusion_times, outlier_lim, sd_lim, win_size, smooth_size, tbw)

% Collects complex cross-spectrum for each epoch, so that statistics can
% be calculated and plotted.
%
% SAMPLE CALL:
% MM_beta_epochs_coh_mtm_collect({'file1.txt','file2.txt'},1000,{'Striatum','Motor
% Ctx.'},7,2,1000,5000,1*2)
%
% INPUTS:
% 'filenames' is a cell of strings, which are the filenames of files 
% containing data for picking beta segments. The data should contain two
% channels, as columns.
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
% 'tbw' is the time bandwidth product.

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

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

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

pd_label = {'pre', 'post'};

for ch = 1:no_channels
    
    for pd = 1:2
        
        fid1 = fopen([file_label, '_', ch_label{ch},...
            '_', pd_label{pd}, '_coh_mtm_', num2str(tbw), 'tbw_r.txt'], 'w');
        
        fid2 = fopen([file_label, '_', ch_label{ch},...
            '_', pd_label{pd}, '_coh_mtm_', num2str(tbw), 'tbw_i.txt'], 'w');
        
        % Loop over filenames.
        
        for fi = 1:length(filenames)
            
            filename = filenames{fi};
            
            beta_listname = [filename,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,'.list'];
            % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
            
            load([beta_listname(1:end-5),'_coh_mtm_', num2str(tbw), 'tbw.mat'])
            
            format = make_format(size(All_coh, 2), 'f');
            
            fprintf(fid1, format, real(All_coh)');
            
            fprintf(fid2, format, imag(All_coh)');
            
        end
        
    end
    
end