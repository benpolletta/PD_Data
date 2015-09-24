function PD_struct = PD_initialize(subject_mat)

load(subject_mat)

PD_struct.folders = folders;

PD_struct.prefixes = prefixes;

PD_struct.basetimes = basetimes;

if exist('infusetimes','var')
    
    PD_struct.infusetimes = infusetimes;

else
    
    PD_struct.infusetimes = nan(size(basetimes));
    
end
    
PD_struct.pd_labels = pd_labels;

PD_struct.chan_labels = chan_labels; PD_struct.no_chans = length(chan_labels);

if exist('striatal_id','var')
    
    PD_struct.striatal_id = striatal_id;

else
    
    PD_struct.striatal_id = [];
    
end

PD_struct.subj_prefix = subject_mat(1:(end - length('_subjects.mat')));

PD_struct.no_folders = length(folders);

try
    
    load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')
    
    PD_struct.sampling_freq = sampling_freq;
    
catch
    
    PD_struct.sampling_freq = 500;
    
end

PD_struct.freqs = 1:200;

PD_struct.bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; PD_struct.no_bands = size(PD_struct.bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(PD_struct.no_bands, 1));

for b = 1:PD_struct.no_bands
   
    band_indices{b} = PD_struct.freqs >= PD_struct.bands(b, 1) & PD_struct.freqs <= PD_struct.bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', PD_struct.bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', PD_struct.bands(b, :));
    
end

PD_struct.band_indices = band_indices;

PD_struct.short_band_labels = short_band_labels;

PD_struct.band_labels = band_labels;

PD_struct.norms = norms; 

PD_struct.no_norms = length(norms);

PD_struct.long_norms = {'', ', Increase Over Baseline Power', ', % Total Power', ', Increase in % Total Power Over Baseline'};

PD_struct.high_type = {'', '_cum'}; 

PD_struct.no_types = length(PD_struct.high_type);

PD_struct.pd_labels = pd_labels; PD_struct.no_pds = length(pd_labels);

if PD_struct.no_pds == 2
    
    PD_struct.long_pd_labels = {'Pre-Infusion', 'Post-Infusion'};
    
elseif PD_struct.no_pds == 4
   
    PD_struct.long_pd_labels = {'Pre-Laser', 'Laser On', '0-5s Post-Laser', '5-10s Post-Laser'};
    
end