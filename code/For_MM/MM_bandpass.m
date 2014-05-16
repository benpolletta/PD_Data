function MM_bandpass(filenames, sampling_freq, conv_length, marker_times)

% 'filenames' is a list of filenames containing data for plotting beta power.
% 'sampling_freq' is the sampling frequency of the data.
% 'conv_length' is the timescale of smoothing of power (will be multiplied
% by sampling freq., so enter in seconds, say, for sampling freq. in Hz).
% 'marker_times' is a list of times to be marked on the resulting power
% time series plots, again to be multiplied by sampling freq. Must have the
% same # of rows as # files; each column contains data for a separate
% sequence of markers.

bands = [10 20; 20 30];

band_names = {'Low Beta','High Beta'};

no_bands = length(band_names);

no_files = length(filenames);

if ~isempty(marker_times)
    
    [r,c] = size(marker_times);
    
    if r~=no_files
        
        if c==no_files
            
            marker_times = marker_times';
            
        else
            
            display('Number of rows of "marker_times" must be the same as the number of files.')
            
            return
            
        end
        
    end
    
    no_markers = size(marker_times,2);
    
else
    
    no_markers = 0;
    
end

for file_no = 1:no_files
    
    filename = filenames{file_no};

    data = load(filename);
    
    if isstruct(data)
        
        fields = fieldnames(data);
        
        data = getfield(data,fields{1});
        
    end
    
    [r,c] = size(data);
    
    if r < c
        
        data = data';
        
    end
    
    [data_length, no_channels] = size(data);
    
    t = (1:data_length)/sampling_freq;
    
    H = nan(data_length,no_channels,no_bands);
    A = nan(data_length,no_channels,no_bands);
    P = nan(data_length,no_channels,no_bands);
    
    if isempty(dir([filename,'_betaHAP.mat']))
        
        clear BP
            
        for b = 1:no_bands
            
            BP = eegfilt(data',sampling_freq,bands(b,1),bands(b,2));
            
            H(:,:,b) = hilbert(BP');
            A(:,:,b) = abs(H(:,:,b));
            P(:,:,b) = angle(H(:,:,b));
            
        end
        
        save([filename,'_betaHAP.mat'],'H','A','P','bands','band_names')
        
    else
        
        load([filename,'_betaHAP.mat'])
        
    end
    
    figure;
    
    [r,c] = subplot_size(no_bands);

    clear A_smooth A_pct_baseline A_plot
    
    flip_length = min(conv_length*sampling_freq,data_length);
    
    for b = 1:no_bands
        
        subplot(c,r,b)
        
        for ch = 1:no_channels
        
            A_flipped = [flipud(A(1:flip_length,ch,b)); A(:,ch,b); flipud(A((end-flip_length+1):end,ch,b))];
            
            A_conv = conv(A_flipped,ones(flip_length,1)/(flip_length),'same');
        
            A_smooth(:,ch) = A_conv((flip_length+1):(end-flip_length));
             
        end
        
        A_plot = A_smooth;
        
        ax = plotyy(t',A_plot(:,1),t',A_plot(:,2));
        
        axis(ax,'tight')
        
        legend({'Striatum','Motor Ctx.'})
        
        hold on
        
        box off

        if ~isempty(marker_times)
            
            marker_time = marker_times(file_no,:);
            
            for m = 1:no_markers
                
                plot([marker_time(m) marker_time(m)]',[min(min(A_plot)) max(max(A_plot))]','r')
                
            end
            
        end
        
        ylabel(sprintf('%s (%g - %g Hz) power',band_names{b},bands(b,1),bands(b,2)))

        xlabel('Time (s)')
        
        if b==1
            
            title(filename)
            
        end
        
    end

    save_as_pdf(gcf,[filename,'_all_channel_data_dec_betaHAP'])

end