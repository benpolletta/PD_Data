datafiles=dir('*.mat');

no_files = length(datafiles);
PD_data = [];

sampling_rate = 20*10^3;

last_index = 0;

for f=1:no_files
    
    load(datafiles(f).name)
    
    if isstruct(data)
        
        data = data.ch1;
        
    end
    
    data = [fliplr(data(1:5*sampling_rate)) data fliplr(data(1:5*sampling_rate))];
    
    data = eegfilt(data, sampling_rate, 0, 300);
    
    data = data((5*sampling_rate + 1):(end - 5*sampling_rate));
    
    data = detrend(data, 'linear');
    
    PD_data(last_index+(1:length(data))) = data;
    
    clear data
    
    last_index = length(PD_data);
    
end

save('All_data.mat','PD_data')