function [mean_E, std_E, normalizer] = PD_powerfreq(subject_mat)

low_freq_lim = 3;
high_freq_lim = 72;

load(subject_mat);
pd_label = {'pre','post'};

normalizer(70, length(folders), length(chan_labels))=0; % Same as setting normalizer = zeros(70, length(folders, length(chan_labels));

for pd=1:length(pd_label)
    
    avg_E=[];
    % pre_E=[];

    for fo = 1:length(folders)-1
        
        folder = folders{fo};
        prefix = prefixes {fo};
        
        outputname=[subject_mat(1:end-12), prefix, '_', pd_label{pd}, '_all_channel_spec_HT', '.mat'];
        
        load(outputname)
        
        for ch = 1:length(chan_labels)
            % pre_E(:,fo,ch)=mean(zscore(log(energy(a:b, :, ch)), 0, 1), 2);
            
            if pd == 1
                
                normalizer(:, fo, ch) = mean(energy(low_freq_lim:high_freq_lim, 1:60*5*fo, ch), 2);
            
            end
            
            pre_E = zeros(70, length(energy), length(chan_labels));
            
            % pre_E(:,:,ch)=energy(a:b,:,ch) ./ repmat(normalizer(:,fo,ch),[1 length(energy) 1]);
            pre_E(:, :, ch) = (energy(low_freq_lim:high_freq_lim, :, ch) - repmat(normalizer(:, fo, ch), [1 length(energy) 1])) ...
                ./ (energy(low_freq_lim:high_freq_lim, :, ch) + repmat(normalizer(:, fo, ch), [1 length(energy) 1]));
            % avg_E(:,ch)=reshape(pre_E(:,fo,ch),b-a+1,1);
            
            avg_E(:, fo, ch) = mean(pre_E(:,:,ch), 2);
                
            figure(fo)
            
            subplot(1, 2, ch)
            
            if pd==1
            
                plot(avg_E(:,fo,ch))
            
            else
                
                plot(avg_E(:,fo,ch),'r')
            
            end
            
            hold on
            
        end

    end
    
    mean_E(:, :, pd) = mean(avg_E, 2);
    std_E(:, :, pd) = std(avg_E, 1, 2) / sqrt(size(avg_E, 2));
   
end
 
 for ch = 1 : length(chan_labels)
     
        figure(fo + 1)
        subplot(1, 2, ch)
        boundedline(1 : length(mean_E), reshape(mean_E(:, ch, :), [70 2]), std_E(:, ch, :))

        title(chan_labels(ch))
        % ylabel(sprintf('log (Power)\n[z-score]'))
        ylabel(sprintf('Power\n(normalized to pre-infusion)'))
        xlabel('Frequency (Hz)')
        
        if ch==2
            legend('Pre-Infusion','Post-Infusion')
        end
 end

 figdir='C:\Users\Administrator\SkyDrive\Boston University\Han Lab\Results PPTs\Carb Paper Figures\';
 save_as_pdf(gcf,[figdir,subject_mat(1:end-4),'_powerfreq'])
