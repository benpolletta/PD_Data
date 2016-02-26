% classify_peaks accepts 4 arguments: a samples x 1 dimensional matrix of
% LFP data, sampling_freq in Hz, prominence (minimum prominence to use for
% peaks) and basetime, which is the time in seconds at which the carbachol
% infusion begins. Clusters is the number of clusters to be used. This
% function returns peak data, which is a matrix
function [peak_indicator, Peak_data] = classify_peaks(folder, prefix, channel_multipliers, sampling_freq, min_prominences, min_secs_apart, basetime, clusters, normalized)
    
load([folder, '/', prefix, '_all_channel_data_dec.mat'])

if isscalar(min_secs_apart)
    min_secs_apart = min_secs_apart*[1 1];
end

if isscalar(clusters)
    clusters = clusters*[1 1];
end

if isscalar(min_prominences)
    min_prominences = min_prominences*[1 1];
end

if isempty(clusters)
    clusters = [2 2];
end

% specify the minimum prominence
if isempty(min_prominences)
    min_prominences = 0;
end

% specify minimum separation in time between the spikes
if isempty(min_secs_apart)
    min_secs_apart = 0.5;
end

if isempty(normalized)
    normalized = 0;
end

if isempty(basetime)
    basetime = 0;
end

rad = floor(min_secs_apart*sampling_freq/2); % 50;

for ch = 1:size(PD_dec, 2)
    
    X = channel_multipliers(ch)*PD_dec(:, ch);
    
    % initialize time array
    t = (1:length(X))/sampling_freq - basetime;
    
    min_samples_apart = min_secs_apart(ch)*sampling_freq;
    
    fig_name = sprintf('%s/%s_chan%d_%dclusters_%.1fHz_spikes_%dprom', folder,...
        prefix, ch, clusters(ch), sampling_freq/min_samples_apart, min_prominences(ch));
    
    if normalized == 1
        fig_name = [fig_name '_normalized'];
    else
        fig_name = [fig_name '_unnormalized'];
    end
    %         find the peaks and their locations within the data
    [peaks, locs] = findpeaks(X, 'MinPeakDistance', min_samples_apart);
    
    %         get peak_details. This function was modified from Ben's original version.
    [peak_widths, peak_prominences, peak_forms] = peak_details(X, locs, min_samples_apart, sampling_freq, normalized, rad);
    
    peak_prom_std = zscore(peak_prominences);
    
    peak_forms(peak_prom_std < min_prominences(ch), :) = [];
    
    peak_forms = peak_forms - repmat(mean(peak_forms, 2), 1, size(peak_forms, 2));
    
    opts = statset('Display','final');
    %
    %     [coeff, X] = pca(peak_forms);
    %
    %     peak_forms_simplified = X*(coeff(:,1:4));
    
    [IDX, C] = kmeans(peak_forms, clusters(ch), 'Replicates', 10, 'Options', opts);
    
    plot_centroids(C, fig_name);
    save_as_pdf(gcf, sprintf('%s_%03d_centroids', fig_name, 0));
    % close
    
    plot_clust_v_time(IDX, locs(peak_prom_std >= min_prominences(ch)), sampling_freq, 20);
    save_as_pdf(gcf, sprintf('%s_%03d_time_histogram', fig_name, 0));
    % close
    
    if rad < sampling_freq
        
        [~, ~, long_peak_forms] = peak_details(X, locs, sampling_freq, sampling_freq, normalized, sampling_freq);
        
        plot_peak_triggered_wt(fig_name, IDX, long_peak_forms, 1:200, linspace(3, 21, length(1:200)), sampling_freq)
        
    else
        
        plot_peak_triggered_wt(fig_name, IDX, peak_forms, 1:200, linspace(3, 21, length(1:200)), sampling_freq)
        
    end
    
    Peak_data = [peaks locs peak_widths peak_prom_std peak_prominences]; %IDX];
    
    Peak_data(peak_prom_std < min_prominences(ch), :) = []; Peak_data = [Peak_data IDX];
    
    plot_peaks_2(fig_name, t, sampling_freq, X, Peak_data, [1 length(X)] + [-min_samples_apart min_samples_apart])
    
    artifact_ids = input('Which cluster number(s) is/are artifacts?');
    
    if ~isempty(artifact_ids)
        
        peak_indicator = scarlet_letter(X, Peak_data, artifact_ids);
        save([prefix '_chan' num2str(ch) '_artifacts.mat'],'peak_indicator','Peak_data');
        
    end
    
    close('all')
    
end

end

% takes these Peak data and ids to remove as arguments. Returns a single
% binary array with 1s at each of these peaks
function outputarray = scarlet_letter(X, Peak_data, artifact_ids)

    outputarray = zeros(size(X));
    
    artifact_locs = [];
    
    for i=1:length(artifact_ids)
        
        curr_locs = Peak_data(Peak_data(:,6) == artifact_ids(i), 2);
        artifact_locs = [artifact_locs; curr_locs];
        
    end
    
    % artifact_matrix = repmat(-rad:1:(rad-1)',length(offset),1);
    % offset_matrix = repmat(offset,1,size(artifact_matrix,2));
    % artifact_matrix = artifact_matrix + offset_matrix;
    % artifact_matrix = sort(unique(artifact_matrix(:)));
    % artifact_matrix(artifact_matrix <= 0) = [];
    % artifact_matrix(artifact_matrix > length(outputarray)) = [];
    % outputarray(artifact_matrix) = 1;
    
    outputarray(artifact_locs) = 1;
    
end

function plot_clust_v_time(IDX, locs, sampling_freq, binlength)
    clusters = unique(IDX);
    clusters = clusters(~isnan(clusters));
    num_clusters = length(clusters);
    xdim = 2;
    ydim = ceil(num_clusters/2);
    figure
    for k=1:num_clusters
        subplot(ydim, xdim,k);
        hist(locs(IDX == k)/sampling_freq,binlength:binlength:(locs(end)/sampling_freq));
        axis tight
        xlabel('Time [s]');
        ylabel('Number of peaks');
        title(['Peak vs. time histogram for cluster ',num2str(k)]);
    end
end


% copy pasta from PD_spikes. edited (see comments)
function [peak_widths, peak_prominences, peak_forms] = peak_details(data, locs, min_distance, sf, normalized, rad)
    % specify sampling frequency
    sampling_freq = sf;
    data_length = length(data);
    % number of peaks is equal to the number of peak locations
    no_peaks = length(locs);
    % initialize peak width and peak prominence arrays
    [peak_widths, peak_prominences] = deal(nan(size(locs)));
    % peak forms is new. This is to store the amplitudes for the samples
    % surrounding each peak
    [peak_forms] = nan(no_peaks,rad*2);

    for p = 1:no_peaks
    %     get location of current peak
        loc = locs(p);
    %   get beginning and end of window that surrounds the peak
        window_start = max(loc - ceil(min_distance/2), 1);

        window_end = min(loc + ceil(min_distance/2), data_length);

        window_indices = window_start:window_end;

    %     get the difference between the peak height and the heights of the
    %     surrounding sample
        dropoff = data(loc) - data(window_indices);
    %     this will mark the beginning of the "peak" period
        peak_start = max(ceil(loc-rad), 1);

    %     this will mark the end of the "peak" period
        peak_end = min(floor(loc+rad), length(data));

    %     this will store the vector 
        peak_vector = data(peak_start:(peak_end-1));

    %     subtract out mean
    %     peak_vector = peak_vector - mean(peak_vector);

    % need stronger normalization. divide through by rms
    %     peak_vector = peak_vector / norm(peak_vector);

    % rms doesn't work. try feature scaling
        if normalized == 1
            disp('Normalizing...');
            peak_vector = (peak_vector-min(peak_vector)) / (max(peak_vector)-min(peak_vector));
        end

    %     add nans if the peak is very close to the beginning or end
        if (loc-rad < 1)
            peak_vector = [nan(abs(loc-rad),1); peak_vector];
        end

        if (length(peak_vector) < rad*2)
            peak_vector = [peak_vector; nan(rad*2-length(peak_vector),1)];
        end

        peak_forms(p,:) = peak_vector';

        peak_prominences(p) = sum(abs(dropoff));

        d_dropoff = diff(dropoff);

        left_endpoint = find(sign(d_dropoff(window_indices < loc)) < 0, 1, 'last');

        offset = length(window_indices)-1-length(d_dropoff(window_indices(1:(end - 1)) > loc));

        right_endpoint = find(sign(d_dropoff(window_indices(1:(end - 1)) > loc)) > 0, 1)+offset;

        try
            peak_widths(p) = right_endpoint - left_endpoint + 1;
            assert(peak_widths(p) >= 0, ['Peak width: ',peak_widths(p)]);
        catch err
            peak_widths(p) = NaN;
            err.stack
        end

    end

end

function plot_peaks_2(fig_name, t, sampling_freq, X, Peak_data, plot_bounds)

    epoch_secs = 20;

    epoch_length = epoch_secs*sampling_freq;

    colors = [0 0 0; 1 0 0; 0 0 1; 0 0.75 0];

    no_epochs = ceil((diff(plot_bounds) + 1)/epoch_length);

    for e = 1:no_epochs

        epoch_start = max(plot_bounds(1) + (e - 1)*epoch_length, 1);

        epoch_end = min(plot_bounds(1) + e*epoch_length - 1, length(X));

        epoch_t = t(epoch_start:epoch_end);

        epoch_data = X(epoch_start:epoch_end, :);

        plot_peaks = nan(length(epoch_t), 1);

        epoch_loc_indices = Peak_data(:, 2) >= epoch_start & Peak_data(:, 2) <= epoch_end;

        epoch_locs = Peak_data(epoch_loc_indices, 2) - epoch_start + 1;

        epoch_loc_t = epoch_t(epoch_locs);

        epoch_peaks = Peak_data(epoch_loc_indices, 1);

        epoch_prominences = Peak_data(epoch_loc_indices, 5);

        epoch_prom_std = Peak_data(epoch_loc_indices, 4);

        plot_peaks(epoch_locs) = epoch_peaks;
        cluster = Peak_data(epoch_loc_indices, 6);

        if mod(e, 3) == 1
            if e < 10
                save_as_pdf(gcf,sprintf('%s_%03d', fig_name, e));
            else    
                save_as_pdf(gcf,sprintf('%s_%03d', fig_name, e));
            end
            % close
            figure
        end

        subplot(3, 1, mod(e, 3) + 3*(mod(e, 3) == 0))

    %     be careful! Added next line...gets rid of second channel's data
        epoch_data = epoch_data(:,1);
        plot(epoch_t, epoch_data)

        hold on

        plot(epoch_t, plot_peaks, 'v')

        uniques = unique(cluster);

        uniques(isnan(uniques)) = [];
        for ind = 1:length(uniques)

            clust = uniques(ind);

            curr_indices = find(cluster == clust);

            text(epoch_loc_t(curr_indices), epoch_peaks(curr_indices)+.05, cellstr(num2str(round(10*clust)/10)),...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Color', colors(mod(clust,length(colors))+1, :), 'FontSize', 12, 'FontWeight', 'bold')
        end

        axis tight

        title(strrep(fig_name(1:10), '_', ' '));

    end
end


function plot_kdist(X,minPts, mouse, ch)
figure
minPts = minPts;
D = squareform(pdist(X));
D = sort(D);
disp(D(1,1));
k_diffs = [];
for i=1:(length(D)-minPts)
    k_diffs(end+1) = D((i+minPts),i);
end
k_diffs = fliplr(sort(k_diffs));
plot(k_diffs);
hold on
xlabel('Index');
ylabel('Distance');
title(['K-distance for mouse ',num2str(mouse),' ch',num2str(ch)]);
save_as_pdf(gcf, ['Kdist_',num2str(mouse),...
    '_ch',num2str(ch)]);
hold off
close

end

function plot_centroids(C, fig_name)
dim = size(C,1);
xdim = 2;
ydim = ceil(dim/2);
figure
yl = [min(C(:)), max(C(:))];

for i=1:dim
    subplot(ydim,xdim,i)
    x = -(length(C(i,:))/2):(length(C(i,:))/2-1);
    plot(x,C(i,:));
    xlim([min(x) max(x)]);
    ylim(yl)
    title(['Cluster number ',num2str(i)]);
    xlabel('Index');
    ylabel('Normalized amplitude');
end

end

function plot_peak_triggered_wt(fig_name, IDX, peak_forms, freqs, no_cycles, sampling_freq)

no_clusters = max(IDX);

data_size = size(peak_forms, 2);

no_freqs = length(freqs);

cycle_lengths = no_cycles.*(sampling_freq./freqs);

wavelets = dftfilt3(freqs, no_cycles, sampling_freq, 'winsize', max(sampling_freq, max(cycle_lengths)));
            
peak_loc = zeros(data_size, 1); peak_loc(round(data_size/2)) = 1;

[Spec_nans, ~] = indicator_to_nans(peak_loc, sampling_freq, freqs, no_cycles, []); Spec_nans = Spec_nans';

segment_length = sampling_freq;

avg_wt_power = nan(no_freqs, data_size, no_clusters);

avg_wt_phase = nan(no_freqs, data_size, no_clusters);

for c = 1:no_clusters
    
    cluster_indices = find(IDX == c);
    
    no_peaks = length(cluster_indices);
    
    wt_power = nan(no_freqs, data_size, no_peaks);
    
    wt_phase = nan(no_freqs, data_size, no_peaks);
    
    parfor p = 1:no_peaks
        
        Spec = nan(no_freqs, data_size);
        
        peak_form = peak_forms(cluster_indices(p), :);
        
        peak_form_reflected = [fliplr(peak_form(1:segment_length)), peak_form, fliplr(peak_form((end - segment_length + 1):end))];
        
        for f = 1:no_freqs
            
            conv_prod = nanconv(peak_form_reflected, wavelets(f,:), 'same');
            
            Spec(f, :) = conv_prod((segment_length + 1):(end - segment_length));
            
        end
        
        % Spec(logical(Spec_nans)) = nan;
        
        wt_power(:, :, p) = abs(Spec);
        
        wt_phase(:, :, p) = angle(Spec);
        
    end
    
    avg_wt_power(:, :, c) = nanzscore(nanmean(wt_power, 3)')';
    
    avg_wt_phase(:, :, c) = nanmean(wt_phase, 3);
    
end

pow_cmax = all_dimensions(@nanmax, avg_wt_power);
pow_cmin = all_dimensions(@nanmin, avg_wt_power);

t = 1000*(((1:data_size) - round(data_size/2))/sampling_freq);

[rows, cols] = subplot_size(no_clusters);

figure

for c = 1:no_clusters
    
    subplot(rows, cols, c)
    
    imagesc(t, freqs, avg_wt_power(:, :, c))
    
    axis xy
    
    % caxis([pow_cmin, pow_cmax])
    
    title(['Cluster # ', num2str(c), ', WT Power'])
    
    xlabel('Time (ms)')
    
    ylabel('Freq. (Hz)')
    
    hold on
    
    plot(t', ones(data_size, 2)*diag([8 30]), 'w:')
    
end

save_as_pdf(gcf, [fig_name, '_wt_pow'])

figure

for c = 1:no_clusters
    
    subplot(rows, cols, c)
    
    for_plot = avg_wt_power(:, :, c);
    
    for_plot(logical(Spec_nans)) = nan;
    
    imagesc(t, freqs, for_plot)
    
    axis xy
    
    % caxis([-pi, pi])
    
    title(['Cluster # ', num2str(c), ', WT Power'])
    
    xlabel('Time (ms)')
    
    ylabel('Freq. (Hz)')
    
end

save_as_pdf(gcf, [fig_name, '_wt_power_nopeak'])

figure

for c = 1:no_clusters
    
    subplot(rows, cols, c)
    
    imagesc(t, freqs, avg_wt_phase(:, :, c))
    
    axis xy
    
    % caxis([pow_cmin, pow_cmax])
    
    title(['Cluster # ', num2str(c), ', WT Phase'])
    
    xlabel('Time (ms)')
    
    ylabel('Freq. (Hz)')
    
end

save_as_pdf(gcf, [fig_name, '_wt_phase'])

end
