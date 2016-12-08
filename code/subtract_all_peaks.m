function [data_subtracted, data] = subtract_all_peaks(subject_mat)

load(subject_mat)

no_subjects = length(folders);

% [r, c] = subplot_size(no_subjects);

present_dir = pwd;

all_fig = figure;

for s = 1:no_subjects
    
    cd (folders{s})
    
    load([prefixes{s}, '_all_channel_data_dec.mat'])
    
    data{s} = PD_dec;
    
    data_subtracted{s} = subtract_peaks(prefixes{s}, '_kmeans', outlier_lims(s), 2);
    
    fig = figure;
    
    plotyy((1:size(PD_dec, 1))', [data{s} data_subtracted{s}], (1:size(PD_dec, 1))', data{s} - data_subtracted{s})
    
    title(folders{s})
    
    save_as_pdf(fig, [prefixes{s}, '_data_subtracted'])
    
    close(fig)
    
    cd (present_dir)
    
    figure(all_fig)
    
    subplot(no_subjects, 1, s)
    
    plotyy((1:size(PD_dec, 1))', [data{s} data_subtracted{s}], (1:size(PD_dec, 1))', data{s} - data_subtracted{s})
    
    title(folders{s})
    
end

save_as_pdf(all_fig, [subject_mat(1:(end - length('_subjects.mat'))), '_data_subtracted.mat'])

close('all')