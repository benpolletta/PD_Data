function PD_list(prefix)

master_fid = fopen([prefix,'_min_master.list'],'w');
epochs_fid = fopen([prefix,'_min_epochs.list'],'w');

for l = -5:14
    
    listname = [prefix,'_min',num2str(l),'.list'];
    
    list_fid = fopen(listname,'w');
    fprintf(master_fid,'%s\n',listname);
    
    start_epoch = l*10;
    
    for e = start_epoch:(start_epoch+9)
   
        epoch_name = [prefix,'_epoch',num2str(e),'.txt'];
        
        fprintf(list_fid,'%s\n',epoch_name);
        fprintf(epochs_fid,'%d\t%d\t%s\n',l,e,epoch_name)
        
    end
       
    fclose(list_fid)
    
end

fclose('all')

master_fid = fopen([prefix,'_5min_master.list'],'w');
epochs_fid = fopen([prefix,'_5min_epochs.list'],'w');

for l = -1:2
    
    listname = [prefix,'_5min',num2str(l),'.list'];
    
    list_fid = fopen(listname,'w');
    fprintf(master_fid,'%s\n',listname);
    
    start_epoch = l*50;
    
    for e = start_epoch:(start_epoch+49)
   
        epoch_name = [prefix,'_epoch',num2str(e),'.txt'];
        
        fprintf(list_fid,'%s\n',epoch_name);
        fprintf(epochs_fid,'%d\t%d\t%s\n',l,e,epoch_name)
        
    end
       
    fclose(list_fid)
    
end

fclose('all')