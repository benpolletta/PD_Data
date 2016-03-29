function Spike_indicator = peak_loader(subj_name, peak_suffix, data_length)

Spike_indicator = zeros(data_length, 2);

if isempty(peak_suffix)
    
    if ~isempty(dir([subj_name, '_peaks.mat']))
        
        load([subj_name, '_peaks.mat'])
        
    elseif ~isempty(dir([subj_name, '_chan1_artifacts.mat'])) || ~isempty(dir([subj_name, '_chan2_artifacts.mat']))
        
        for ch = 1:2
            
            if ~isempty(dir([subj_name, '_chan', num2str(ch), '_artifacts.mat']))
                
                load([subj_name, '_chan', num2str(ch), '_artifacts.mat'])
                
                Spike_indicator(:, ch) = peak_indicator;
            
        
            end
        
        end
        
    end

elseif strcmp(peak_suffix, '_kmeans')
        
    if ~isempty(dir([subj_name, '_chan1_artifacts.mat'])) || ~isempty(dir([subj_name, '_chan2_artifacts.mat']))
        
        for ch = 1:2
            
            if ~isempty(dir([subj_name, '_chan', num2str(ch), '_artifacts.mat']))
                
                load([subj_name, '_chan', num2str(ch), '_artifacts.mat'])
                
                Spike_indicator(:, ch) = peak_indicator;
                
            elseif ~isempty(dir([subj_name, '_peaks.mat']))
                
                peaks = load([subj_name, '_peaks.mat']);
                
                peaks_Si = peaks.Spike_indicator;
                
                Spike_indicator(:, ch) = peaks_Si(:, ch);
                
            end
            
        end
   
    end
end