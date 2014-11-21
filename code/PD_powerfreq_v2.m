function [mean_E,std_E,normalizer] = PD_powerfreq_v2(subject_mat)
sampling_freq=1000;
a=3;
b=72;

load(subject_mat);
pd_label = {'pre','post'};

normalizer(70,length(folders),length(chan_labels))=0;

for pd=1:length(pd_label)
    
    avg_E=[];
    pre_E=[];
    
    for fo = 1:length(folders)-1
        
        folder = folders{fo};
        prefix = prefixes {fo};
        
        outputname=[pwd,'/',subject_mat(1:end-12),prefix,'_',pd_label{pd},'_all_channel_spec_HT', '.mat'];
        
        load(outputname)
        
        for ch=1
            pre_E=zeros(70,length(energy),length(chan_labels));
            pre_E(:,:,ch)=energy(a:b,:,ch);
            
            beta_listname = [folder,'/',prefix,'_ch1_beta_post_7out_2sd_333.3333win_5000smooth_win.list'];
            
            betapower = [folder,'/',prefix,'_all_channel_data_dec_HAP.mat'];
            
            % Loading block numbers, epoch start and end indices.
            [blocks, beta_starts, beta_ends, epochs, epoch_starts, epoch_ends] = text_read(beta_listname,'%f%f%f%f%f%f%*[^\n]');
            
            %             find(max(epochs)
            betatime = beta_starts(find(epochs==max(epochs))-max(epochs)+1);
            
        end
        
        for ch=1:length(chan_labels)
            %             pre_E(:,fo,ch)=mean(zscore(log(energy(a:b,:,ch)),0,1),2);
            if pd==1
                normalizer(:,fo,ch) = mean(energy(a:b,1:30*sampling_freq,ch),2);
%                 normalizer(:,fo,ch) = mean(energy(a:b,:,ch),2);
            end
            
            pre_E=zeros(70,length(energy),length(chan_labels));
            
            %             pre_E(:,:,ch)=energy(a:b,:,ch) ./ repmat(normalizer(:,fo,ch),[1 length(energy) 1]);
            pre_E(:,:,ch)=(energy(a:b,:,ch) - repmat(normalizer(:,fo,ch),[1 length(energy) 1]))./ (energy(a:b,:,ch) + repmat(normalizer(:,fo,ch),[1 length(energy) 1]));
            %             avg_E(:,ch)=reshape(pre_E(:,fo,ch),b-a+1,1);
            
            
            figure(fo)
            subplot(1,2,ch)
            if pd==1
                avg_E(:,fo,ch) = mean(pre_E(:,1:30*sampling_freq,ch),2);
                plot(avg_E(:,fo,ch))
                hold on
                
                figure(fo+10)
                subplot(2,2,ch)
                imagesc(1:size(pre_E(:,1:30*sampling_freq,ch),2),1:size(pre_E(:,1:30*sampling_freq,ch),1),pre_E(:,1:30*sampling_freq,ch))
                axis xy
            else
%                 avg_E(:,fo,ch)=mean(pre_E(:,betatime:betatime+30*sampling_freq,ch),2);
                avg_E(:,fo,ch)=mean(pre_E(:,infusetimes(fo)*sampling_freq:(infusetimes(fo)+30)*sampling_freq,ch),2);
                plot(avg_E(:,fo,ch),'r')
                hold on
                
                figure(fo+10)
                subplot(2,2,ch+2)
                imagesc(1:size(pre_E(:,infusetimes(fo)*sampling_freq:(infusetimes(fo)+30)*sampling_freq,ch),2),1:size(pre_E(:,infusetimes(fo)*sampling_freq:(infusetimes(fo)+30)*sampling_freq,ch),1),pre_E(:,infusetimes(fo)*sampling_freq:(infusetimes(fo)+30)*sampling_freq,ch))
                axis xy
            end
            hold on
            
            
        end
        
    end
    
    mean_E(:,:,pd)=mean(avg_E,2);
    std_E(:,:,pd)=std(avg_E,1,2)/sqrt(size(avg_E,2));
    
    
end


for ch=1:length(chan_labels)
        figure(fo+1)
        subplot(1,2,ch)
        boundedline(1:length(mean_E),reshape(mean_E(:,ch,:),[70 2]),std_E(:,ch,:))

        title(chan_labels(ch))
%         ylabel(sprintf('log (Power)\n[z-score]'))
        ylabel(sprintf('Power\n(normalized to pre-infusion)'))
        xlabel('Frequency (Hz)')
        
        if ch==2
            legend('Pre-Infusion','Post-Infusion')
        end
 end

 figdir='C:\Users\Administrator\SkyDrive\Boston University\Han Lab\Results PPTs\Carb Paper Figures\';
 save_as_pdf(gcf,[figdir,subject_mat(1:end-4),'_powerfreq_v2'])