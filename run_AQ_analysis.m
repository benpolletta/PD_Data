present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

for fo = 1:length(folders)
    
    folder = folders{fo};
   
    cd (folder)
    
    PD_concatenate
    
    files = dir('*.mat');
    
    for fi = 1:length(files)
        
        filename = files(fi).name;
        
        load(filename)
        
        if ~isstruct(data)
            
            data = struct('ch1',data);
           
        end
        
        filenames{fi} = [filename(1:end-4), '_HT.mat'];
       
        spec_v2(pwd, filename, 1)
        
    end
    
    parkpop(pwd, filenames, 1)
    
    cd (present_dir)
    
end