function BetaTimes = run_AQ_analysis

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

basetimes = [300 1200 1800 600 1800];

infusetimes = [390 240 300 450 510];

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
        
        save_as_pdf(gcf,[prefix,'_PPBP.fig'])
        
    end
    
    cd (present_dir)
    
end