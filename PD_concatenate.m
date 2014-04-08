datafiles=dir('*.mat');

no_files=length(datafiles);
PD_data=[];

last_index=0;

for f=1:no_files
    
    load(datafiles(f).name)
    
    PD_data(last_index+(1:length(data)))=data;
    
    last_index=length(PD_data);
    
end

save('PD_data.mat','PD_data')