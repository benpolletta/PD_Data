function delta_theta_ratio_plot

load('PD_delta_theta_ratio')

load('M1_groups')

M1_minus_index = ~M1_indices(:, 1);
M1_plus_index = ~M1_indices(:, 2);

M1_minus_dt = delta_theta(M1_minus_index, 2);
M1_plus_dt = delta_theta(M1_plus_index, 2);

figure

plot(cumsum(ones(6, 2), 2), [M1_plus_dt M1_minus_dt], 'o')

xlim([0 3])

set(gca, 'XTick', [1 2], 'XTickLabel', {'M1+', 'M1-'})

ylabel('\delta/\theta Ratio')

box off

set(gca, 'FontSize', 16)

[~, p] = ttest2(M1_plus_dt, M1_minus_dt);
    
sigstar([1 2], p)
    
end