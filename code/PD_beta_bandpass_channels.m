function PD_beta_bandpass_channels

present_dir = pwd;

% load('initial_subjects.mat')

load('st_m1_subjects.mat')

bands = [10 20; 20 30];
band_names = {'Low Beta','High Beta'};
no_bands = length(band_names);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    infusetime = infusetimes(fo);

    load([folder,'/',prefix,'_all_channel_data_dec.mat'])
    
    t = (1:length(PD_dec))/sampling_freq;
    
    H = nan(size(PD_dec,1),2,no_bands);
    A = nan(size(PD_dec,1),2,no_bands);
    P = nan(size(PD_dec,1),2,no_bands);
    
    if isempty(dir([folder,'/',prefix,'_all_channel_data_dec_betaHAP.mat']))
        
        BP = nan(length(PD_dec),2);
            
        for b = 1:no_bands
            
            BP = eegfilt(PD_dec',sampling_freq,bands(b,1),bands(b,2));
            
            H(:,:,b) = hilbert(BP');
            A(:,:,b) = abs(H(:,:,b));
            P(:,:,b) = angle(H(:,:,b));
            
        end
        
        save([folder,'/',prefix,'_all_channel_data_dec_betaHAP.mat'],'H','A','P','bands','band_names')
        
    else
        
        load([folder,'/',prefix,'_all_channel_data_dec_betaHAP.mat'])
        
    end
    
    figure;
    
    [r,c] = subplot_size(no_bands);

    clear A_smooth A_pct_baseline A_plot
    
    for b = 1:no_bands
        
        subplot(c,r,b)
        
        for ch = 1:2
        
            A_flipped = [flipud(A(1:10*sampling_freq,ch,b)); A(:,ch,b); flipud(A((end-10*sampling_freq+1):end,ch,b))];
            
            A_conv = conv(A_flipped,ones(20*sampling_freq,1)/(20*sampling_freq),'same');
        
            A_smooth(:,ch) = A_conv((10*sampling_freq+1):(end-10*sampling_freq));
             
        end
            
        A_pct_baseline = 100*A_smooth*diag(1./mean(A_smooth(t<basetime,:))) - 100;
        
        A_plot = A_smooth;%A_pct_baseline;
        
        plot(t',A_plot)
        
        legend({'Striatum','Motor Ctx.'})
        
        hold on
        
        axis tight
        
        box off
        
        plot([basetime basetime]',[min(min(A_plot)) max(max(A_plot))]','r')
        plot([basetime+infusetime basetime+infusetime]',[min(min(A_plot)) max(max(A_plot))]','r')

        ylabel(sprintf('%s (%g - %g Hz) power',band_names{b},bands(b,1),bands(b,2)))
%         ylabel({[band_names{b},' power'];'percent change'})

        xlabel('Time (s)')
        
        if b==1
            
            title(folder)
            
        end
        
    end

    save_as_pdf(gcf,[folder,'/',prefix,'_all_channel_data_dec_betaHAP'])
%     save_as_pdf(gcf,[folder,'/',prefix,'_all_channel_data_dec_betaHAP_pct'])

end