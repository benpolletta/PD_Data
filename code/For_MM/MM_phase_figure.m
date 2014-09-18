function MM_phase_figure(filenames, ~, chan_labels, ~, outlier_lim, sd_lim, win_size, smooth_size, tbw)

% Combines phase plots from roseplot and coherence analyses, for final
% figure for the paper.
%
% SAMPLE CALL: MM_phase_figure({'file1.txt','file2.txt'},1000,{'Striatum','Motor
% Ctx.'},7,2,1000,5000,1*2)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

% Creating name to call the analysis by - either name of single file in
% 'filenames', or 'All'. NB: will be over-written each time you run a
% multi-file analysis, so should be renamed after running.

if length(filenames) == 1
    
    file_label = filenames{1};
    
else

    file_label = 'All';
    
end

f_bins = 9.5:1:30.5; no_f_bins = length(f_bins) - 1;

f_pairs = nchoosek(1:no_f_bins, 2); no_f_pairs = size(f_pairs, 1);

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

f_center_indices = f_centers <= 25;

% c_order = [linspace(1,0,no_f_bins); abs(linspace(1,0,no_f_bins)-.5); linspace(0,1,no_f_bins)]';

f_labels = textscan(num2str(f_centers), '%s', 'delimiter', ' ');
f_labels = cellstr(f_labels{1});
f_labels = f_labels(1:2:end);

f = 1000*(0:(win_size - 1))/win_size;

f_indices = f <= 26 & f >= 9;

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

figure;

rad_deg = 180/pi;

%% Plotting phase angle by frequency, pre and post.

for ch = 1:4
    
    for ch1 = 1:2
        
        figure
        
        load([file_label, '_', par_name, '_', ch_label{ch}, '_by_', ch_label{ch1}, '_beta_ri_rose_dp_group.mat']);
        
        subplot(2, 1, 1)
        
        conf_mat = reshape(conf_mat, size(conf_mat, 1), 1, size(conf_mat, 2));
        
        conf_mat = repmat(conf_mat, [1 2 1]);
        
        h = boundedline(f_centers(f_center_indices)', rad_deg*angle(MR_mat(f_center_indices, :)), rad_deg*conf_mat(f_center_indices, :, :));
        
        set(h, 'Marker', 's')
        
        axis tight
        
        hold on
        
        plot(f(f_indices)', zeros(length(f(f_indices)), 1), '--k')
            
        legend(period_label, 'Location', 'SouthEast')
        
        axis tight
        
        xlabel('Frequency (Hz)')
        
        title({[file_label, ' High Beta Blocks']; ['Mean Phase Angle (', chan_labels{1}, ' - ', chan_labels{2}, ') by ', chan_labels{1}, ' Freq.']})
        
        %% Plotting coherence.
        
        subplot(2, 1, 2)
        
        coh_listname = ['All_', ch_label{ch}, '_post_coh_mtm_', num2str(tbw), 'tbw_phase.mat'];
        
        load(coh_listname)
        
        h = boundedline(f(f_indices)', rad_deg*mean_data(f_indices, :), rad_deg*std_data(f_indices, :, :));
        
        set(h, 'Marker', 's')
        
        hold on
        
        plot(f(f_indices)', zeros(length(f(f_indices)), 1), '--k')
        
        legend(period_label)
        
        axis tight
        
        xlabel('Frequency (Hz)')
        
        title({[chan_labels{ch}, ' High Beta Blocks'];'Phase of Coherence'})
        
    end
    
end

save_as_pdf(gcf, [file_label, '_phase_figure'])

end