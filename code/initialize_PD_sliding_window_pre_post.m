%% General stuff

load('seven_bands')

if ~exist('band_index', 'var'), band_index = 4; end

data_labels_struct = init_data_labels(freqs, no_cycles, bands, 'data_field', 'data_subtracted', 'band_index', band_index)

subjects_mat_cell = {'st_m1_subjects.mat', 'st_m1_ali_subjects.mat', 'st_m1_ali2_subjects.mat'};

chan_labels = {'Striatum', 'M1'};

filename = 'STR_w_M1';

if ~exist('window_length', 'var'), window_length = 150; end

sliding_window_cell{1} = window_length*data_labels_struct.sampling_freq{1}*[1 1];

%% Specific stuff.

if ~exist('analysis', 'var'), analysis = 'granger'; end

switch analysis
    
    case 'granger'
        
        % Granger causality analysis:
        
        sliding_window_cell{2} = [2 2];
        
        fcn_handle = @mvgc_analysis;
        
        varargin = {[], '', 1};
        
        shuffle_struct = init_struct({'shuffle_dims', 'no_shuffles', 'combine_pairs', 'varargin'},...
            {1, 1000, @combine_dimension, {2, [1 2]}});
        
    case 'pmtm'
        
        % PMTM:
        
        sliding_window_cell{2} = [1 1];
        
        fcn_handle = @pmtm_detrend;
        
        varargin = {window_length, [], data_labels_struct.sampling_freq{1}};
        
    case 'arspec'
        
        % AR spectrum:
        
        sliding_window_cell{2} = [1 1]; % {[500 500], [1 1]}; %
        
        fcn_handle = @mvgc_analysis;
        
        varargin = {[], '', 3};
        
    case 'arxspec'
        
        % AR cross-spectrum:
        
        sliding_window_cell{2} = [2 2]; % = {[150*500 150*500], [2 2]};
        
        fcn_handle = @mvgc_analysis;
        
        varargin = {[], '', 3};
        
    case 'pac'
        
        sliding_window_cell{2} = [1 1]; % = {[150*500 150*500], [2 2]};
        
        fcn_handle = @PAC;
        
        varargin = {data_labels_struct.sampling_freq{1}, [.25:.25:10 11:30], [15:30 32.5:2.5:80 85:5:200]};
        
        shuffle_struct = init_struct({'shuffle_dims', 'no_shuffles', 'combine_pairs', 'varargin'},...
            {1, 1000, @pass_input, {}});
        
end

function_name = func2str(fcn_handle);