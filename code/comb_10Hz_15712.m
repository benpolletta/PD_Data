function comb_10Hz_15712

load('15712/15712_all_channel_data.mat')

% sampling_freq = 20*10^3;

if isempty(dir('15712/15712_all_channel_data_uncombed.mat'))

    save('15712/15712_all_channel_data_uncombed.mat')

end
    
% Comb filtering at 10.31 Hz & harmonics.

for f = 1:20
    
    notch_freq = 10.31*f;
    
    [n,Wn] = buttord(2*(notch_freq + [-1 1])/sampling_freq, 2*(notch_freq + 1.25*[-1 1])/sampling_freq, 1, 5);
    
    [z, p, k] = butter(n, Wn, 'stop'); [sos, g] = zp2sos(z, p, k); h = dfilt.df2sos(sos, g);
    
    for ch = 1:2
        
        PD_data(:, ch) = filtfilthd(h, PD_data(:, ch), 'reflect');
        
    end
    
end

save('15712/15712_all_channel_data.mat', 'PD_data', 'sampling_rate')