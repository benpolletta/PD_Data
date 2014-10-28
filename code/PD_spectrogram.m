function PD_spectrogram(subject_mat)
sampling_freq=1000;
a=1;
b=100;

load(subject_mat);

% pd_label = {'pre','post'};
% period_label = {'Pre-Infusion','Post-Infusion'};
% 
% pd_label = {'6ohda'};
% period_label = {'6OHDA'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    prefix = prefixes {fo};
    
    base_index = basetimes(fo)*sampling_freq;
    beta_index = (betatimes(fo))*sampling_freq+base_index;
    subj_name = [folder,'/',prefix];
    
    load([pwd,'/',subj_name,'_all_channel_data_dec.mat'])
    
    for pd=1:length(pd_label)
        
        for i=1:length(chan_labels)
            channel=['ch' num2str(i)];
            dummydata=sampling_freq*5; %5sec of dummy data
            if pd == 1
                % tempdata = [PD_dec(dummydata:-1:1,i); PD_dec(1:base_index,i); PD_dec(base_index:-1:base_index-dummydata+1,i)]';
                    tempdata = [PD_dec(dummydata:-1:1,i); PD_dec(:,i); PD_dec(end:-1:end-dummydata+1,i)]';
            else if pd == 2
                    tempdata = [PD_dec(dummydata+base_index-1:-1:base_index,i); PD_dec(base_index+1:beta_index,i); PD_dec(beta_index:-1:beta_index-dummydata+1,i)]';
                end
            end
            
            SegPoints=5000;
            
            % Generates a spectrogram of the recording - Figure (24).
            fprintf(['Generating a HT-Spectrogram for ',prefix,' ',channel,'\n'])
            phase=[];
            energy=[];
            
            % figure(fo)
            [phase,energy] = basic_HT_improved_x11(tempdata, sampling_freq, SegPoints, dummydata, a, b);
            title([num2str(folder),' - ',chan_labels{i},' ',period_label{pd}])
            % save_as_pdf(gcf,[subject_mat(1:end-12),prefix,'_',pd_label{pd},'_',channel,'_spec'])
            
            
            PD_spec.phase(:,:,i) = single(phase);
            PD_spec.energy(:,:,i) = single(energy);
            PD_spec.LFP(:,i) = tempdata(dummydata:end-dummydata-1);
            clear phase energy
        end
        
        
        outputname = [subject_mat(1:end-12), prefix, '_', pd_label{pd}, '_all_channel_spec_HT', '.mat'];
        save (outputname, '-struct','PD_spec','-v6')
        clear PD_spec
        
    end
    
end