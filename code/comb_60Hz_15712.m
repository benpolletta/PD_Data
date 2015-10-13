function comb_60Hz_15712

load('15712/15712_all_channel_data.mat')

% sampling_freq = 20*10^3;

if isempty(dir('15712/15712_all_channel_data_un-60Hz-combed.mat'))

    save('15712/15712_all_channel_data_un-60Hz-combed.mat')

end

% Comb filtering at 60 Hz & harmonics.

for f = 1:3
    
    notch_freq = 60*f;
    
    [n,Wn] = buttord(2*(notch_freq + 2*[-1 1])/sampling_rate, 2*(notch_freq + 2.5*[-1 1])/sampling_rate, 1, 5);
    
    [z, p, k] = butter(n, Wn, 'stop'); [sos, g] = zp2sos(z, p, k); h = dfilt.df2sos(sos, g);
    
    for ch = 1:2
        
        PD_data(:, ch) = filtfilthd(h, PD_data(:, ch), 'reflect');
        
    end
    
end

save('15712/15712_all_channel_data.mat', 'PD_data', 'sampling_rate')