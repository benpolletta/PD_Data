function PD_phase_diff

sampling_freq = 1000;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

basetimes = [300 1200 1800 600 1800];

infusetimes = [390 240 300 450 510];

bands = [1 4; 4 12; 10 30; 30 60; 60 90; 90 110; 120 180];
band_names = {'delta','theta','beta','lgamma','hgamma','HFO'};
% bands = [10 20; 20 30];
% band_names = {'Low Beta','High Beta'};
no_bands = length(band_names);

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    infusetime = infusetimes(fo);
    
    if isempty(dir([folder,'/',prefix,'_all_channel_data_dec_P_diff.mat']))
        
        load([folder,'/',prefix,'_all_channel_data_dec_HAP.mat'])
        
        P_diff = reshape(diff(unwrap(P),[],2),size(P,1),no_bands);
        
        P_diff = angle(exp(sqrt(-1)*P_diff));
        
        zero_indicator = P_diff > -pi/6 & P_diff < pi/6;
        
        P_diff_smooth = nan(size(P_diff));
        
        zi_smooth = nan(size(zero_indicator));
        
        zi_peaks = cell(no_bands,1);
        
        zero_indices = cell(no_bands,1);
        
        for b = 1:no_bands
            
            P_diff_flipped = [flipud(P_diff(1:1*sampling_freq,b)); P_diff(:,b); flipud(P_diff((end-1*sampling_freq+1):end,b))];
            
            P_diff_conv = conv(exp(sqrt(-1)*P_diff_flipped),ones(.5*sampling_freq,1)/(.5*sampling_freq),'same');
            
            P_diff_smooth(:,b) = angle(P_diff_conv((1*sampling_freq+1):(end-1*sampling_freq)))/pi;
            
            zi_smooth(:,b) = conv(single(zero_indicator(:,b)),ones(.5*sampling_freq,1)/(.5*sampling_freq),'same');
            
            [~,zi_peak_indicator] = spaced_peaks(zi_smooth(:,b),.5*sampling_freq,0);
            
            zero_indices{b} = zi_peak_indicator & (zi_smooth(:,b) >= .5);
            
        end
        
        save([folder,'/',prefix,'_all_channel_data_dec_P_diff.mat'],'P_diff','P_diff_smooth','zero_indicator','zi_smooth','zero_indices')
        
    else
        
        load([folder,'/',prefix,'_all_channel_data_dec_P_diff.mat'])
        
    end
    
    t = (1:size(P_diff,1))/sampling_freq;
    
    figure;
    
    [r,c] = subplot_size(no_bands);
        
    for b = 1:no_bands
        
        subplot(c,r,b)
        
        plot(t',P_diff_smooth(:,b))
        
        hold on
        
        axis tight
        
        box off
        
        plot(t',zi_smooth(:,b),'g')
        
        plot(t(zero_indices{b})',zeros(sum(zero_indices{b}),1),'k*')
        
        y_lims = [min(min(P_diff_smooth(:,b))) max(max(P_diff_smooth(:,b)))]';
        
        plot([basetime basetime]',y_lims,'r')
        plot([basetime+infusetime basetime+infusetime]',y_lims,'r')

        ylabel(sprintf('%s (%g - %g Hz) Phase Difference (pi)',band_names{b},bands(b,1),bands(b,2)))

        xlabel('Time (s)')
        
        if b==1
            
            title(folder)
            
        end
        
    end

    save_as_pdf(gcf,[folder,'/',prefix,'_all_channel_data_dec_P_diff'])

end