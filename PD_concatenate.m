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
    
    PD_data(last_index+(1:length(data))) = data;
    
    clear data
    
    last_index = length(PD_data);
    
end

save('All_data.mat','PD_data','sampling_rate')