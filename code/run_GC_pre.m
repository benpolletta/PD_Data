function run_GC_pre

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

basetimes = [300 1200 1800 600 1800];

infusetimes = [390 240 300 450 510];

load('BetaTimes')

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    PD_concatenate_channels(prefix,1)
    
    PD_decimate_channels(prefix,1)

    cd (present_dir)
    
end

get_beta