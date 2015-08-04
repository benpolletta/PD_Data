function [X_vals, Y_vals] = nans_to_patch(times, freqs, nan_indicator)

no_freqs = length(freqs);

no_patches = sum(diff([0; nan_indicator(:, 1); 0]) == 1);

max_patch_length = 2*no_freqs;

[left_border_indices, right_border_indices] = deal(cell(no_freqs));

for f = 1:no_freqs
    
    left_border_indices{f} = find(diff([0; nan_indicator(:, f); 0]) == 1);
    
    right_border_indices{f} = find(diff([0; nan_indicator(:, f); 0]) == -1);
    
    if f >= 2 && sum(left_border_indices{f}) > sum(left_border_indices{f - 1})
    
        max_patch_length = max_patch_length + 2*(no_freqs - f + 1);
    
    end
    
end

[X_vals, Y_vals] = deal(cell(no_patches, 1));

border_indicator = abs(diff([zeros(1, no_freqs); nan_indicator; zeros(1, no_freqs)]));

for p = 1:no_patches
    
    patch_border = border_indicator(left_border_indices{1}(p):right_border_indices{1}(p), :);
    
    [J_vals, I_vals] = find(patch_border');
    
    % IJ_vals = [I_vals J_vals];
    % 
    % IJ_vals = sortrows(IJ_vals, 1);
    % 
    % I_vals = IJ_vals(:, 1); J_vals = IJ_vals(:, 2);
    
    X_vals{p} = times(I_vals);
    
    Y_vals{p} = freqs(J_vals);
    
end
    
end
    