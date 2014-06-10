function BetaTimes = run_AQ_analysis

present_dir = pwd;

load('subjects.mat')

BetaTimes = nan(length(folders),4);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    datalist = [prefix,'_datafiles.list'];
    
    if isempty(dir(datalist))
    
        file_fid = fopen(datalist,'w');
        
        for digit=1:9
            
            file = dir([prefix,'*',num2str(digit),'.mat']);
            
            if ~isempty(file)
                
                fprintf(file_fid,'%s\n',file.name);
                
            end
            
        end
        
        fclose(file_fid);
    
    end
    
    file_fid = fopen(datalist,'r');
    
    clear filenames
    
    files = textscan(file_fid,'%s');
    files = files{1};
    
    for fi = 1:length(files)
        
        filename = files{fi};
        
        filenames{fi} = [filename(1:end-4), '_HT.mat'];
        
        if isempty(dir(filenames{fi}))
       
            spec_v2(pwd, filename, 1)
    
        end
            
    end
    
    if isempty(dir('ParkPop.mat'))
        
        parkpop(pwd, filenames, 1)
        
    else
        
        chData = load('ParkPop.mat');
        
        BetaTimes(fo,:) = Park_pop_betapower_v3(chData.ch1, basetimes(fo), infusetimes(fo));
        
        save_as_pdf(gcf,[prefix,'_PPBP'])
        
    end
    
    cd (present_dir)
    
end