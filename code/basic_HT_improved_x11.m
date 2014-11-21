function [phase,energy] = basic_HT_improved_x11(LFP, srate, SegPoints, dummypoints, low_freq_lim, high_freq_lim, plot_opt)

% Made by Niccolo Talei Franzesi, the decimation proceedure done by Mike Henninger.
% You will not get stuck in while loop as long as (6*safetymargin/maxfilt<1
% and safetymargin^2<3*maxfilt)

% maxfiltlength=4000;
bandwidth = 3;
% n = last channel/trial of interest

% filteredsignal=LFP;
    
Segment_No = floor(length(LFP)/SegPoints);

no_freqs = high_freq_lim - low_freq_lim + 1; 

freqs = nan(no_freqs, 1);
phase = nan(no_freqs, Segment_No*SegPoints);
energy = nan(no_freqs, Segment_No*SegPoints);

freq_index = 1;

for lowfreq = low_freq_lim:high_freq_lim
    
    highfreq = lowfreq + bandwidth;
    
    freqs(freq_index) = mean(lowfreq, highfreq);
    
    freq_index = freq_index + 1;
    
    for j = 1:Segment_No
        
        % if j == Segment_No
        %   tempanaLFP=LFP((j-1)*SegPoints+1:end);
        % else
        
        Segment_start = max(1, (j - 1)*SegPoints + 1 - dummypoints);
        Segment_end = min(length(LFP), j*SegPoints + dummypoints);
        
        temp_LFP = detrend(LFP(Segment_start:Segment_end)); % Detrends (brings to zero mean) the filtered signal by subtracting a linear best fit (not just average!).
        
        % end
        
        temptempsignal = eegfilt(temp_LFP, srate, lowfreq, highfreq);
        
        tempsignal = temptempsignal(dummypoints + 1:end - dummypoints); %% To get rid of edge artifact 0.1s beginning and end.
        
        filt_det = detrend(tempsignal); % Detrends (brings to zero mean) the filtered signal by subtracting a linear best fit (not just average!).
        
        phase(lowfreq + 2, (j-1)*SegPoints + 1:length(tempsignal)) = angle(hilbert(filt_det));  % If we didn't downsample then we just spit it out.
        energy(lowfreq + 2, (j-1)*SegPoints + 1:length(tempsignal)) = abs(hilbert(filt_det));
        
        clear tempsignal;
        clear filt_det;
        warning off all;
        
    end
    
end

if plot_opt == 1
    
    figure
    
    t = length(energy)/srate;
    
    % climLower=min(median(energy));
    % climUpper=max(median(energy));
    % clims = [climLower/2 climUpper*2];
    % imagesc(energy, clims)
    
    imagesc(t, freqs, zscore(energy')')
    
    axis xy
    
    % axis ([0 length(energy)*0.001 5 200]);
    
    xlabel ('Time (s)');
    
    ylabel ('Frequency (Hz)');
    
end



