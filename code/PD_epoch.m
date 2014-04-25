function PD_epoch(prefix, injection_sec)

load([prefix,'_all_data_dec.mat'])

all_mean = mean(PD_dec);
all_std = std(PD_dec);

t_secs = 1:length(PD_dec)/sampling_freq;

epoch_length = 6;

no_pre = ceil(injection_sec/epoch_length);
no_post = ceil(t_secs(end)/epoch_length);

for e = -no_pre:1:(no_post-1)
   
    start_sec = injection_sec - (e+1)*epoch_length;
    
    end_sec = injection_sec - e*epoch_length;
    
    start_index = max(start_sec*sampling_freq + 1, 1);
    
    end_index = min(end_sec*sampling_freq, length(PD_dec));
    
    data = PD_dec(start_index:end_index);
    
    fid = fopen([prefix,'_epoch',num2str(e),'.txt'],'w');
    
    fprintf(fid,'%f\n',data);
    
    fclose(fid);
    
end