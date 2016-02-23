% classify_peaks accepts 4 arguments: a samples x 1 dimensional matrix of
% LFP data, sampling_freq in Hz, prominence (minimum prominence to use for
% peaks) and basetime, which is the time in seconds at which the carbachol
% infusion begins. Clusters is the number of clusters to be used. This
% function returns peak data, which is a matrix
function Peak_data = classify_peaks_mike(X, sampling_freq, min_prominences, min_secs_apart, basetime, clusters, normalized)

    if isempty(clusters)
        clusters = 2;
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

    rad = min_secs_apart*sampling_freq/2; % 50;
    
    % initialize time array
    t = (1:length(X))/sampling_freq - basetime;

    min_samples_apart = min_secs_apart*sampling_freq;

    fig_name = sprintf('%dclusters%.1fHz_spikes_%dprom',clusters, ...
            sampling_freq/min_samples_apart, min_prominences);
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
    
    peak_forms(peak_prom_std < min_prominences, :) = [];
    
    peak_forms = peak_forms - repmat(mean(peak_forms, 2), 1, size(peak_forms, 2));

    opts = statset('Display', 'final');

    [IDX, C] = kmeans(peak_forms, clusters, 'Replicates', 5, 'Options', opts);

    plot_centroids(C);

    plot_clust_v_time(IDX, locs(peak_prom_std >= min_prominences), sampling_freq, 20);

    Peak_data = [peaks locs peak_widths peak_prom_std peak_prominences]; % IDX];

    Peak_data(peak_prom_std < min_prominences, :) = []; Peak_data = [Peak_data IDX];

    plot_clust_v_time(Peak_data(:, end), Peak_data(:, 2), sampling_freq, 20);

    plot_peaks_2(fig_name, t, sampling_freq, X, Peak_data, [1 length(X)] + [-min_samples_apart min_samples_apart])

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
        peak_start = max(floor(loc-rad),1);

    %     this will mark the end of the "peak" period
        peak_end = min(floor(loc+rad),length(data));

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

            text(epoch_loc_t(curr_indices), epoch_peaks(curr_indices)+.2, cellstr(num2str(round(10*clust)/10)),...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Color', colors(mod(clust,length(colors))+1, :), 'FontSize', 12, 'FontWeight', 'bold')
        end

        axis tight

        title(fig_name(1:6))
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

function plot_centroids(C)
dim = size(C,1);
xdim = 2;
ydim = ceil(dim/2);
figure
for i=1:dim
    subplot(ydim,xdim,i)
    x = -(length(C(i,:))/2):(length(C(i,:))/2-1);
    plot(x,C(i,:));
    xlim([min(x) max(x)]);
    title(['Cluster number ',num2str(i)]);
    xlabel('Index');
    ylabel('Normalized amplitude');
end

end
