function PD_plot_artifacts(prefix)

outliers = text_read([prefix,'_outliers.list'],'%s%*[^\n]');

no_outliers = length(outliers);

[r,c] = subplot_size(no_outliers);

if no_outliers > 0
    
    figure;
    
    for o = 1:no_outliers
        
        subplot(r,c,o)
        
        data = load(outliers{o});
        
        t = (1:length(data))/1000;
        
        plot(t,data,'k')
        
        title(outliers{o})
        
    end
    
    save_as_pdf(gcf,[prefix,'_outliers']);
    
end