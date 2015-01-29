function collect_striatal_freq(norm)

subject_matnames = {'st_m1', 'st_stn'};

channels = [1 2];

load('st_m1_subjects.mat', 'pd_labels')

no_pds = length(pd_labels);

fid = nan(no_pds, 1);

for pd = 1:no_pds
    
    fid(pd) = fopen(['STR_beta_block_freqs', norm, '_ch1_', pd_labels{pd}, '.txt'], 'w');
    
end

for pd = 1:no_pds
    
    for s = 1:2
        
        freq_data = load([subject_matnames{s}, '_beta_block_freqs', norm, '_ch', num2str(channels(s)), '_', pd_labels{pd}, '.txt']);
        
        format = make_format(size(freq_data, 2), 'f');
        
        fprintf(fid(pd), format, freq_data');
        
    end
    
end

fclose('all')