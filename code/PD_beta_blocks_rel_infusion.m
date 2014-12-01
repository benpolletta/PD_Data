function PD_beta_blocks_rel_infusion(subject_mat, sd_lim)

% if isscalar(win_size)
% 
%     par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];
% 
% elseif length(win_size) == 2
%    
%     par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size(1)), 'to', num2str(win_size(2)),'win_',num2str(smooth_size),'smooth'];
%     
% else
%     
%     display('win_size must be a scalar (lower limit of beta epoch length) or an interval (lower and upper limits).')
%     
%     return
%     
% end
    
load(subject_mat)

sampling_freq = 1000;

norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);

[BP_high_dps, BP_high_cum_dps] = deal(nan(length(folders), 6, no_norms, 2, 2));

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    base_index = basetimes(fo)*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    for n = 1:no_norms
        
        BP_data = load([subj_name,'_wt_BP.mat'], ['BP', norms{n}]);
        
        BP_data = getfield(BP_data, ['BP', norms{n}]);
        
        t = (1:size(BP_data, 1))/sampling_freq;
        
        BP_high = nan(size(BP_data));
        
        pd_limits = [1 min(length(t), base_index); min(length(t), base_index + 1) min(length(t), base_index + 1500*sampling_freq)];
        
        pd_lengths = diff(pd_limits, [], 2) + 1;
        
        for ch = 1:2
            
            for b = 1:size(BP_high, 2)
                
                high_cutoff = mean(BP_data(1:pd_limits(1, 2), b, ch)) + sd_lim*std(BP_data(1:pd_limits(1, 2), b, ch));
                
                BP_high(:, b, ch) = BP_data(:, b, ch) >= high_cutoff;
                
            end
            
        end
        
        BP_high_cum = BP_high;
        
        BP_high_cum(cumsum(BP_high, 2) > 1) = 0;
        
        for pd = 1:size(pd_limits,1)
            
            for ch = 1:2
                
                BP_high_dps(fo, :, ch, pd, n) = sum(BP_high(pd_limits(pd,1):pd_limits(pd,2), :, ch))/pd_lengths(pd);
                    
                BP_high_cum_dps(fo, :, ch, pd, n) = sum(BP_high_cum(pd_limits(pd,1):pd_limits(pd,2), :, ch))/pd_lengths(pd);
                
            end
            
        end
        
        save([subj_name, '_', num2str(sd_lim), 'sd_BP', norms{n}, '_high.mat'], 'BP_high', 'BP_high_cum')
        
    end
          
end

save([subject_mat(1:(end - length('_subjects.mat'))), '_', num2str(sd_lim), 'sd_BP_high_dps.mat'], 'BP_high_dps', 'BP_high_cum_dps')