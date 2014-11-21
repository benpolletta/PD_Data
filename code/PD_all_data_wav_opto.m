function PD_all_data_wav_opto(subjects_mat)

subjects_struct = load(subjects_mat);

folders = subjects_struct.folders;

prefixes = subjects_struct.prefixes;

freqs = 1:200; no_freqs = length(freqs);

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];

no_cycles = linspace(3, 21, no_freqs);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    wav_inner_opto(folder, prefix, freqs, no_cycles, bands)
    
end