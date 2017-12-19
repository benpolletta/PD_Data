function compare_motor_groups_PAC(norm_struct, significance)

close('all')

%% Looping over normalizations.

if strcmp(norm_struct, 'all')

    % Plotting normalization by baseline or by shuffles.

    whos = {'baseline', 'shuffle'};

    hows = {'', 'subtract', 'zscore'};

    for w = 1:length(whos)

        for h = 1:length(hows)

            norm_struct = struct('who', whos{w}, 'how', hows{h});

            compare_motor_groups_PAC(norm_struct)

        end

    end

    % Plotting sequential normalization by baseline and by shuffles.

    whos = {'shuffle_baseline', 'baseline_shuffle'};

    hows = {'', 'subtract', 'zscore'};

    for w = 1:length(whos)

        for h = 1:length(hows)

            for h1 = 1:length(hows)

                norm_struct = struct('who', whos{w});

                norm_struct.how = {hows{h}, hows{h1}};

                compare_motor_groups_PAC(norm_struct)

            end

        end

    end

    return

end

%% Plotting PAC with xPlt.

clear SW_xPlt

analysis = 'pac'; window_length = 10;

initialize_PD_sliding_window_pre_post

window_time_cell = cellfun(@(x,y) x/y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

pd_names = {'pre', 'post'}; no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods

    pd_label = [pd_label, '_', pd_names{period}];

end

if nargin == 0, norm_struct = []; end
if isempty(norm_struct)
    norm_struct = struct('who', 'shuffle_baseline');
    norm_struct.how = {'zscore', ''};
end

norm_label = ['_', norm_struct.who];

if isfield(norm_struct, 'how')

    if ~isempty(norm_struct.how)

        if ischar(norm_struct.how)

            norm_label = [norm_label, '_', norm_struct.how];

        elseif iscellstr(norm_struct.how)

            if any(cellfun(@(x) ~isempty(x), norm_struct.how))

                for c = 1:length(norm_struct.how), norm_label = [norm_label, '_', norm_struct.how{c}]; end

            end

        end

    end

end

load('M1_groups')

M1_increased_index{2} = M1_increased_index{2} & All_index{2};

M1_not_increased_index{2} = M1_not_increased_index{2} & All_index{2};

groups_plotted = {M1_increased_index, M1_not_increased_index}; % {All_index};

for group = 1:length(groups_plotted)

    group_name = [make_sliding_window_analysis_name([filename, groups_plotted{group}{1}, pd_label,...
        '_band', num2str(data_labels_struct.band_index), make_label('bands', legnth(data_labels_struct.bands), 7, 'back')], function_name,...
        window_time_cell, 2, varargin{:}), norm_label];

    SW_group = load([group_name, '_recordingspacked.mat']);

    SW_group = squeeze(SW_group.SW_RecordingsPacked);

    group_axis = nDDictAxis;
    group_axis.name = 'Group';
    group_axis.values = {groups_plotted{group}{1}(2:end)};

    SW_group.axis(end + 1) = group_axis;

    if exist('SW_xPlt', 'var')

        SW_xPlt = SW_xPlt.merge(SW_group);

    else

        SW_xPlt = SW_group;

    end

end

SW_xPlt.getaxisinfo

function_handles = {@xp_tight_subplot_adaptive, @xp_compare_3D};
dimensions = {{'Channel', 'Period'}, {'Group'}};
function_arguments = {{},{@ranksum, significance}};
recursivePlot(SW_xPlt, function_handles, dimensions, function_arguments)

function_arguments = {{},{[], significance}};
recursivePlot(SW_xPlt, function_handles, dimensions, function_arguments)

SW_post_minus_pre = norm_axis_by_value(SW_xPlt, 'Period', 'pre', 'subtract');
SW_post_minus_pre = SW_post_minus_pre.axissubset('Period', 'post');

function_arguments = {{},{@ranksum, significance}};
dimensions = {{'Channel'}, {'Group'}};
recursivePlot(SW_post_minus_pre, function_handles, dimensions, function_arguments)

function_arguments = {{},{[], significance}};
recursivePlot(SW_post_minus_pre, function_handles, dimensions, function_arguments)

save_all_figs(['compare_motor_groups_PAC', norm_label, '_p', num2str(significance)])

% frequencies = SW_RecordingsPacked.meta.matrix_dim_1.values;
% freq_index = frequencies <= freq_limit;
% SW_RecordingsPacked.data = cellfun(@(x) x(freq_index, :), SW_RecordingsPacked.data, 'UniformOutput', false);
% SW_RecordingsPacked.meta.matrix_dim_1.values = frequencies(freq_index);
%
% SW_cross = SW_RecordingsPacked.axissubset('Channel 1', channel_labels{1});
% SW_cross = SW_cross.axissubset('Channel 2', channel_labels{2});
%
% ax(4, group) = subplot(no_measures_plotted, 2, 2*3 + group);
%
% xp_comparison_plot_2D(squeeze(SW_cross), str2func(test_flag), 2*p_val, [], 1)
%
% legend off
%
% if group == 1
%
%     ylabel('Cross-Spectrum', 'FontSize', 10)
%
% end
%
% sync_axes(ax(4, :))
