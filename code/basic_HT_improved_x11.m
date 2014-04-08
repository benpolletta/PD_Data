function [phase,energy] = basic_HT_improved_x11(LFP,srate,SegPoints, dummypoints, a,b,filename)
%%

%made by Niccolo Talei Franzesi, the decimation proceedure done by Mike Henninger.
%you will not get stuck in while loop as long as (6*safetymargin/maxfilt<1
%and safetymargin^2<3*maxfilt)



%maxfiltlength=4000;
bandwidth=3;
%n = last channel/trial of interest

%note: in labeling, need to make sure I know if it's channel 1,2 or-i guess
%that would be obvious-3

%    filteredsignal=LFP;
    for lowfreq=a:1:b;
        %lowfreq
        highfreq=lowfreq+bandwidth;
        
        SegmentN=ceil(length(LFP)/SegPoints);
       
        for j=1:SegmentN-2
           % if j==SegmentN
           % tempanaLFP=LFP((j-1)*SegPoints+1:end);
           % else                
            tempanaLFP=LFP((j-1)*SegPoints+1:j*SegPoints+2*dummypoints);
           % end
            temptempsignal = eegfilt(tempanaLFP, srate, lowfreq,highfreq);
            tempsignal=temptempsignal(dummypoints+1:end-dummypoints); %%to get rid of artifact 0.1s beginning and end.
            filt_det = detrend(tempsignal); % detrends (brings to zero mean) the filtered signal by subtracting a linear best fit (not just average!)
            phase(lowfreq+2,(j-1)*SegPoints+1:(j-1)*SegPoints+length(tempsignal)) = angle(hilbert(filt_det));  %if we didn't downsample then we just spit it out
            energy(lowfreq+2,(j-1)*SegPoints+1:(j-1)*SegPoints+length(tempsignal)) = abs(hilbert(filt_det));
            clear tempsignal;
            clear filt_det;
            warning off all;
        end
    end
%     figure (23)
%     %skip the zero lines
%     climLower=min(median(energy));
%     climUpper=max(median(energy));
%     clims = [climLower/2 climUpper*2];
%     imagesc(energy, clims)
%     axis xy
%     fname = filename;
%     title(fname)
%     %axis([xmin xmax ymin ymax])
%     %xlabel('time 9.1ms')
%     ylabel('frequency (Hz)')
    
    figure (24)
    x=0:0.001:length(energy)*0.001;
    y=0:200;
    climLower=min(median(energy));
    climUpper=max(median(energy));
    clims = [climLower/2 climUpper*2];
    imagesc(energy, clims)
    imagesc(x, y, energy, clims)
    axis xy
    axis ([0 length(energy)*0.001 5 200]);
    xlabel ('time (s)');
    ylabel ('Frequency (Hz)');


%     figure_filename = ['C:\Documents and Settings\Administrator\My Documents\auto import&start analysis\' fname '.jpg'];



