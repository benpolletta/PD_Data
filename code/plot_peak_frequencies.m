function plot_peak_frequencies(data_labels_struct, norm_struct)

%% Initializing.

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
    
load([filename, '_band', num2str(data_labels_struct.band_index), norm_label, '_peak_frequencies.mat'])

%% Creating 'Cartesian product' of measures.

max_value_axis = nDDictAxis; max_value_axis.name = 'Observable'; max_value_axis.values = {'Max. Value'};
freq_axis = nDDictAxis; freq_axis.name = 'Observable'; freq_axis.values = {'Frequency'};

SW_Power.axis(end + 1) = max_value_axis;
SW_Freq.axis(end + 1) = freq_axis;

max_xp = merge(SW_Power, SW_Freq);
max_xp = xp_matrix_transpose(squeeze(max_xp.packDim('Recording')));
max_xp = max_xp.unpackDim(2, length(max_xp.axis) + 1);
max_xp.axis(max_xp.findaxis('Measure')).name = 'Measure 1';
max_xp = max_xp.repmat(max_xp.axis('Measure 1').values, 'Measure 2', 2);
max_xp_copy = max_xp.permute([2 1 3:length(max_xp.axis)]);

max_xp.data = cellfun(@(x, y) [x y], max_xp.data, max_xp_copy.data, 'UniformOutput', false);

%% Plotting by group.

dim_order = {'Measure 1', 'Measure 2', 'Period', 'Band', 'Observable'};

load('M1_groups')

M1_increased_index{2} = M1_increased_index{2}(All_index{2});
M1_not_increased_index{2} = M1_not_increased_index{2}(All_index{2});
All_index{2} = All_index{2}(All_index{2});

groups_plotted = {All_index, M1_increased_index, M1_not_increased_index};

for group = 1:length(groups_plotted)
    
    group_name = groups_plotted{group}{1};
    
    group_max_xp = max_xp;
    group_max_xp.data = cellfun(@(x) x(groups_plotted{group}{2}, :), max_xp.data, 'UniformOutput', false);
    
    close('all')
    
    function_handles = {@xp_tight_subplot_adaptive, @xp_scatter_w_corr};
    dimensions = {1:5,0}; % {{'Measure 1', 'Measure 2', 'Period', 'Band', 'Observable'}, {}}; %
    function_arguments = {{dim_order},{}};
    [~, titles] = recursivePlot_2(group_max_xp, function_handles, dimensions, function_arguments);
    
    for d = 3:5
        titles = strrep(titles, [dim_order{d}, ': '], '');
    end
    titles = strrep(titles, ' ', '_');
    
    for fig = 1:length(titles)
        save_as_pdf(fig, ['scatter_', group_name, titles{fig}(1:(end - 1))])
    end
    
    close('all')
    
    function_handles = {@xp_tight_subplot_adaptive, @xp_compare_1D};
    function_arguments = {{dim_order},{}};
    [~, titles] = recursivePlot_2(group_max_xp, function_handles, dimensions, function_arguments);
    
    for d = 3:5
        titles = strrep(titles, [dim_order{d}, ': '], '');
    end
    titles = strrep(titles, ' ', '_');
    
    for fig = 1:length(titles)
        save_as_pdf(fig, ['compare_', group_name, titles{fig}(1:(end - 1))])
    end
    
end

end