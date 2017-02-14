function data_labels_struct = init_data_labels(freqs, no_cycles, bands, varargin)

[freqs, no_cycles, bands, BP_suffix] = init_freqs(freqs, no_cycles, bands);

keys = {'BP_suffix', 'peak_suffix', 'data_suffix', 'epoch_secs', 'pd_handle', 'data_field', 'band_index', 'sampling_freq'};

values = {BP_suffix, '_kmeans', '_peakless', 150, '_by_STR', 'PD_dec', 4, {500, 1}};

data_labels_struct = init_struct(keys, values, varargin{:});