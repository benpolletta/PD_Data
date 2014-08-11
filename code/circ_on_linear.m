function [line_params, all_params] = circ_on_linear(l, c, spacing, plot_opt)

rotations = linspace(-pi, pi, spacing);

corrs = nan(size(rotations));

slopes = nan(size(rotations));

intercepts = nan(size(rotations));

for r = 1:length(rotations)
    
    c_rotated = angle(exp(sqrt(-1)*(c + rotations(r))));
    
    corrs(r) = corr(l, c_rotated)^2;
   
    line_params = polyfit(l, c_rotated, 1);
    
    slopes(r) = line_params(1); 
    
    intercepts(r) = line_params(2) - rotations(r);
    
end

all_params = [rotations; corrs; slopes; intercepts];

[~, max_index] = max(corrs);

line_params = [slopes(max_index) intercepts(max_index)];

if plot_opt == 1
   
    figure
    
    subplot(1, 3, 1)
    
    plot(rotations', corrs')
    
    axis tight
    
    title('Correlation')
    
    subplot(1, 3, 2)
    
    plot(rotations', slopes')
    
    axis tight
    
    title('Slope')
    
    subplot(1, 3, 3)
    
    plot(rotations', intercepts')
    
    axis tight
    
    title('Intercept')
    
end