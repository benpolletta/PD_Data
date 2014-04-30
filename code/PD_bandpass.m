present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

basetimes = [300 1200 1800 600 1800];

infusetimes = [390 240 300 450 510];

bands = [1 4; 4 12; 15 30; 30 60; 60 90; 90 110; 120 180];
band_names = {'delta','theta','beta','lgamma','hgamma','HFO'};
no_bands = length(band_names);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    infusetime = infusetimes(fo);

    load([folder,'/',prefix,'_all_data_dec.mat'])
    
    t = (1:length(PD_dec))/sampling_freq;
    
    if isempty(dir([folder,'/',prefix,'_all_data_dec_HAP.mat']))
        
        BP = nan(length(PD_dec),no_bands);
        
        for b = 1:no_bands
            
            BP(:,b) = eegfilt(PD_dec',sampling_freq,bands(b,1),bands(b,2));
            
        end
        
        H = hilbert(BP);
        A = abs(H);
        P = angle(H);
        
        save([folder,'/',prefix,'_all_data_dec_HAP.mat'],'H','A','P','bands','band_names')
        
    else
        
        load([folder,'/',prefix,'_all_data_dec_HAP.mat'])
        
    end
    
    figure;
    
    [r,c] = subplot_size(no_bands);

    for b = 1:no_bands
        
        subplot(r,c,b)
        
        A_smooth = conv(A(:,b),ones(sampling_freq,1)/sampling_freq,'same');
        
        A_pct_baseline = 100*A_smooth/mean(A_smooth(t<basetime)) - 100;
        
        A_plot = A_pct_baseline;
        
        plot(t',A_plot)
        
        hold on
        
        axis tight
        
        box off
        
        plot([basetime basetime]',[min(A_plot) max(A_plot)]','g')
        plot([basetime+infusetime basetime+infusetime]',[min(A_plot) max(A_plot)]','r')

        % ylabel([band_names{b},' power'])
        ylabel({[band_names{b},' power'];'percent change'})
        
        if b==1
            
            title(folder)
            
        end
        
    end

%     save_as_pdf(gcf,[folder,'/',prefix,'_all_data_dec_HAP'])
    save_as_pdf(gcf,[folder,'/',prefix,'_all_data_dec_HAP_pct'])

end