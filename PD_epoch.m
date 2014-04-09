load('All_data_dec.mat')

t_secs = 1:length(PD_dec)/sampling_freq;

injection_sec = 300;

t_rel_injection = t_secs - injection_sec;

epoch_length = 6;

for e = -50:150
   
    start_sec = injection_sec + e*epoch_length;
    
    end_sec = injection_sec + (e+1)*epoch_length;
    
    start_index = start_sec*sampling_freq + 1;
    
    end_index = end_sec*sampling_freq;
    
    data = PD_dec(start_index:end_index);
    
    fid = fopen(['epoch',num2str(e),'.txt'],'w');
    
    fprintf(fid,'%f\n',data);
    
    fclose(fid)
    
end