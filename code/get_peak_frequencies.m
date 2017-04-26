function get_peak_frequencies(data_labels_struct, norm_struct)

if nargin < 1, data_labels_struct = []; end
if isempty(data_labels_struct)
    seven_bands_struct = load('seven_bands');
    data_labels_struct = init_data_labels(seven_bands_struct.freqs, seven_bands_struct.no_cycles, seven_bands_struct.bands, 'data_suffix', '');
end

load('seven_bands')

if nargin < 2, norm_struct = []; end
if isempty(norm_struct), norm_struct = struct('who', 'baseline', 'how', ''); end

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how'), if ~isempty(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how]; end, end

filename = 'STR_w_M1';

load('STR_M1_subjects.mat')

pd_names = {'pre', 'post'}; 
pd_label = '';
for pd = 1:length(pd_names)
    pd_label = [pd_label, '_', pd_names{pd}];
end

load('missing_2'), load('M1_groups')

%% Getting frequencies from spectra.

load([filename, pd_label,'_band', num2str(data_labels_struct.band_index),...
    data_labels_struct.BP_suffix, data_labels_struct.peak_suffix, data_labels_struct.pd_handle, '_pct_spectrum'])

frequencies = SW_xPlt.meta.matrix_dim_1.values;
SW_xPlt = squeeze(mean_over_axis(SW_xPlt, 'Window_Dim_1'));

[~, ~, band_labels] = band_max(SW_xPlt.data{1}, frequencies, bands);

[sw_Power, sw_Freq] = deal(SW_xPlt);

[sw_Power.data, sw_Freq.data] = cellfun(@(x) band_max(x, frequencies, bands(1:end, :)), SW_xPlt.data, 'UniformOutput', 0);

SW_Power = sw_Power; SW_Freq = sw_Freq;

SW_Pow_Ch_vals = SW_Power.axis('Channel').values;
for v = 1:length(SW_Pow_Ch_vals)
    new_values{v} = [SW_Pow_Ch_vals{v}, ' Pow.'];
end

SW_Power.axis(SW_Power.findaxis('Channel')).values = new_values;
SW_Power.axis(SW_Power.findaxis('Channel')).name = 'Measure';
SW_Freq.axis(SW_Freq.findaxis('Channel')).values = new_values;
SW_Freq.axis(SW_Freq.findaxis('Channel')).name = 'Measure';

%% Getting frequencies from PLV.

load([filename, pd_label,'_band', num2str(data_labels_struct.band_index),...
    data_labels_struct.BP_suffix, data_labels_struct.peak_suffix, data_labels_struct.pd_handle, '_pct_PLV'])

SW_xPlt = squeeze(mean_over_axis(SW_xPlt, 'Window_Dim_1'));

[sw_Power, sw_Freq] = deal(SW_xPlt);

[sw_Power.data, sw_Freq.data] = cellfun(@(x) band_max(x, frequencies, bands(1:end, :)), SW_xPlt.data, 'UniformOutput', 0);

measure_axis = nDDictAxis;
measure_axis.name = 'Measure';
measure_axis.values = {'PLV'};

sw_Power.axis(end + 1) = measure_axis; sw_Freq.axis(end + 1) = measure_axis;

sw_Power = sw_Power.alignAxes(SW_Power); SW_Power = SW_Power.merge(sw_Power);

sw_Freq = sw_Freq.alignAxes(SW_Freq); SW_Freq = SW_Freq.merge(sw_Freq);

SW_Freq = squeeze(SW_Freq.packDim('Recording'));
SW_Freq.data = cellfun(@(x) x(:, All_index{2}), SW_Freq.data, 'UniformOutput', false);
SW_Freq = SW_Freq.unpackDim(2, [], 'Recording');
SW_Freq.axis(SW_Freq.findaxis('Recording')).values = folders(All_index{2});

SW_Power = squeeze(SW_Power.packDim('Recording'));
SW_Power.data = cellfun(@(x) x(:, All_index{2}), SW_Power.data, 'UniformOutput', false);
SW_Power = SW_Power.unpackDim(2, [], 'Recording');
SW_Power.axis(SW_Power.findaxis('Recording')).values = folders(All_index{2});

%% Getting frequencies from Cross-Spectrum.
    
pd_names = {'pre', 'post'}; no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods
    
    pd_label = [pd_label, '_', pd_names{period}];
    
end

if nargin < 6, norm_struct = []; end
if isempty(norm_struct), norm_struct = struct('who', 'baseline', 'how', ''); end

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how'), if ~isempty(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how]; end, end

channel_labels = {'Striatum', 'Motor Ctx.'};

analysis = 'arxspec'; window_length = 10; band_index = data_labels_struct.band_index;

initialize_PD_sliding_window_pre_post

window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

band_name = [make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index)], function_name,...
    window_time_cell, 2, varargin{:}), norm_label];

load([band_name, '_recordingspacked.mat'])

frequencies = SW_RecordingsPacked.meta.matrix_dim_1.values;

SW_cross = SW_RecordingsPacked.axissubset('Channel 1', channel_labels{1});
SW_cross = SW_cross.axissubset('Channel 2', channel_labels{2});

SW_cross = squeeze(SW_cross.unpackDim(2, 1, [], folders(All_index{2})));
[sw_Power, sw_Freq] = deal(SW_cross);
[sw_Power.data, sw_Freq.data] = cellfun(@(x) band_max(x, frequencies, bands(1:end, :)), SW_cross.data, 'UniformOutput', 0);

measure_axis = nDDictAxis;
measure_axis.name = 'Measure';
measure_axis.values = {'X-Spec.'};

sw_Power.axis(end + 1) = measure_axis; sw_Freq.axis(end + 1) = measure_axis;

sw_Power = sw_Power.alignAxes(SW_Power); SW_Power = SW_Power.merge(sw_Power);

sw_Freq = sw_Freq.alignAxes(SW_Freq); SW_Freq = SW_Freq.merge(sw_Freq);

%% Getting frequencies for Granger Causality.
    
pd_names = {'pre', 'post'}; no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods
    
    pd_label = [pd_label, '_', pd_names{period}];
    
end

if nargin < 6, norm_struct = []; end
if isempty(norm_struct), norm_struct = struct('who', 'baseline', 'how', ''); end

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how'), if ~isempty(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how]; end, end

short_chan_labels = {'Str.', 'M1'};
    
analysis = 'granger'; window_length = 10;

initialize_PD_sliding_window_pre_post

window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

band_name = [make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index)], function_name,...
    window_time_cell, 2, varargin{:}), norm_label];

load([band_name, '_recordingspacked.mat'])

channel_loc = [2 1];

for direction = 1:2
    
    channel_loc = fliplr(channel_loc);
    
    SW_direction = SW_RecordingsPacked.axissubset('Channel From', channel_labels{channel_loc(1)});
    SW_direction = SW_direction.axissubset('Channel To', channel_labels{channel_loc(2)});
    
    SW_direction = squeeze(SW_direction.unpackDim(2, 1, [], folders(All_index{2})));
    [sw_Power, sw_Freq] = deal(SW_direction);
    [sw_Power.data, sw_Freq.data] = cellfun(@(x) band_max(x, frequencies, bands(1:end, :)), SW_direction.data, 'UniformOutput', 0);
    
    measure_axis = nDDictAxis;
    measure_axis.name = 'Measure';
    measure_axis.values = {sprintf('PDC, %s->%s', short_chan_labels{channel_loc})};
    
    sw_Power.axis(end + 1) = measure_axis; sw_Freq.axis(end + 1) = measure_axis;
    
    sw_Power = sw_Power.alignAxes(SW_Power); SW_Power = SW_Power.merge(sw_Power);
    
    sw_Freq = sw_Freq.alignAxes(SW_Freq); SW_Freq = SW_Freq.merge(sw_Freq);
    
end

band_axis = nDDictAxis; band_axis.name = 'Band'; band_axis.values = band_labels;
SW_Power.meta.matrix_dim_1 = band_axis;
SW_Freq.meta.matrix_dim_1 = band_axis;
    
save([filename, '_band', num2str(data_labels_struct.band_index), norm_label, '_peak_frequencies.mat'], 'SW_Freq', 'SW_Power', 'data_labels_struct', 'norm_struct')

end

function [bm, bmfreq, band_labels] = band_max(data, frequencies, bands)
    
no_bands = size(bands, 1);

if isnumeric(frequencies)
    
    [bm, bmfreq] = deal(nan(no_bands, 1)); band_labels = {};
    
    for b = 1:no_bands
        
        band_labels{b} = sprintf('%d-%d Hz', bands(b, :));
        
        band_indicator = frequencies >= bands(b, 1) & frequencies <= bands(b, 2);
        
        [bm(b), bmi] = max(data(band_indicator));
        
        band_freqs = frequencies(band_indicator);
        
        bmfreq(b) = band_freqs(bmi);
        
    end
    
end

end