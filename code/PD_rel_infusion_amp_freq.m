function PD_rel_infusion_amp_freq(subject_mat)

load(subject_mat)

no_channels = length(chan_labels);

period_label = {'Pre-Infusion','Post-Infusion'};

no_periods = length(period_label);

all_amp_name = [subject_mat(1:(end-length('_subject.mat'))),'_amp.txt'];

all_amp_data = load(all_amp_name);

all_pds = all_amp_data(:, 1);

all_freq_name = [subject_mat(1:(end-length('_subject.mat'))),'freq.txt'];

all_freq_data = load(all_freq_name);

a_norm_labels = {'Power','Power (Rank)'};

an_labels = {'a','ar'};

no_a_norms = length(a_norm_labels);

f_norm_labels = {'Freq.','Freq. (Rank)'};

fn_labels = {'f','fr'};

no_f_norms = length(f_norm_labels);

for a_norm = 1:no_a_norms
    
    for f_norm = 1:no_f_norms
        
        figure
        
        for pd = 1:no_periods
            
            for ch = 1:no_channels
                
                for ch1 = 1:no_channels
                    
                    subplot(no_channels*no_channels, no_periods, (ch - 1)*no_channels*no_periods + (ch1 - 1)*no_periods + pd)
                    
                    A_pd = all_amp_data(all_pds == pd, (a_norm - 1)*no_periods + ch + 1);
                    
                    F_pd = all_freq_data(all_pds == pd, (f_norm - 1)*no_periods + ch1 + 1);
                    
                    r = nancorr(A_pd, F_pd);
                    
                    [histogram, centers] = hist3([A_pd F_pd], [50 50]);
                    
                    colormap('default')
                    
                    imagesc(centers{1}, centers{2}, histogram/sum(sum(histogram)))
                    
                    h = colorbar;
                    
                    title({[chan_labels{ch},' by ',chan_labels{ch1},', ',period_label{pd}];['Correlation = ',num2str(r)]})
                    
                    ylabel(h, 'Prop. Observed')
                    
                    axis xy
                    
                    xlabel([chan_labels{ch}, ' ', a_norm_labels{a_norm}])
                    
                    ylabel([chan_labels{ch1}, ' ', f_norm_labels{f_norm}])
                    
                end
                
            end
        
        end

        save_as_pdf(gcf, [subject_mat(1:(end-length('_subject.mat'))),'_',an_labels{a_norm},'_by_',fn_labels{f_norm}])
    
    end
    
end