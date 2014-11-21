function PD_spectrogram(subject_mat)

low_freq_lim = 1;
high_freq_lim = 100;

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
    
    for pd = 1:length(pd_label)
        
        for ch = 1:length(chan_labels)
            
            channel=['ch' num2str(ch)];
            
            dummydata=sampling_freq*5; %5sec of dummy data
            
            if pd == 1
                
                % tempdata = [PD_dec(dummydata:-1:1, ch); PD_dec(1:base_index, ch); PD_dec(base_index:-1:(base_index - dummydata + 1), ch)]';
                tempdata = [PD_dec(dummydata:-1:1, ch); PD_dec(:, ch); PD_dec(end:-1:end-dummydata+1, ch)]';
                
            elseif pd == 2
                
                tempdata = [PD_dec((dummydata + base_index - 1):-1:base_index, ch); PD_dec(base_index+1:beta_index, ch);...
                    PD_dec(beta_index:-1:(beta_index - dummydata + 1), ch)]';
                
            end
            
            SegPoints = 5000;
            
            % Generates a spectrogram of the recording - Figure (24).
            fprintf(['Generating a HT-Spectrogram for ',prefix,' ',channel,'.\n'])
            
            % figure(fo)
            
            [phase,energy] = basic_HT_improved_x11(tempdata, sampling_freq, SegPoints, dummydata, low_freq_lim, high_freq_lim, 0);
            
            % title([num2str(folder),' - ',chan_labels{ch},' ',period_label{pd}])
            % save_as_pdf(gcf,[subject_mat(1:end-12),prefix,'_',pd_label{pd},'_',channel,'_spec'])
            
            
            PD_spec.phase(:, :, ch) = single(phase(:, (dummydata + 1):(end - dummydata - 1)));
            PD_spec.energy(:, :, ch) = single(energy(:, (dummydata + 1):(end - dummydata - 1)));
            PD_spec.LFP(:, ch) = tempdata((dummydata + 1):(end - dummydata - 1));
            
            clear phase energy
        
        end
        
        outputname = [subject_mat(1:end-12), prefix, '_', pd_label{pd}, '_all_channel_spec_HT', '.mat'];
        save (outputname, '-struct', 'PD_spec', '-v7')
        clear PD_spec
        
    end
    
end