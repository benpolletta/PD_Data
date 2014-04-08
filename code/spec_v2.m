function [] = spec_v2(folderName,fileName,channels)
rate=20000;
a=1;
b=200;

% matfile={DataAnalysisDirectory, filenumber, '.mat'};

% [DataDirectory, DataAnalysisDirectory] = SetDirectory_temp(folderName);

% matfile=abf_to_mat(folderName, fileName); %%convert .abf to .mat for wave_clus
% load(matfile);

present_dir = pwd;

cd (folderName)

load(fileName)

for i=1:channels
    
    channel=['ch' num2str(i)];
    if isstruct(data)
        data = data.(channel);
    end
    
    dummydata=rate*5; % 5 sec. of dummy data.
    %Zeropadding=zeros(1,dummydata);
    %tempdata=horzcat(Zeropadding, data, Zeropadding);
    tempdata = horzcat(data(dummydata:-1:1), data, data(end:-1:end-dummydata+1)); % Dealing with edge effects by reflection.
    
    filteredsignal = eegfilt(tempdata,rate,0,300);  % This performs a lowpass at the safe cutoff frequency.
    clear tempdata;
    
    %Hz60signal = eegfilt(data,rate,59.5,60.5);
    %LFP = downsample(filteredsignal-Hz60signal,20);
    %clear Hz60signal;
    
    downfactor=rate/1000;
    srate = 1000;
    LFP = downsample(filteredsignal,downfactor); % Downsample by a rate of 20;
    clear filteredsignal;
    
    dummypoints=dummydata/downfactor; % Number of "dummy" (or reflected) data points.
    SegPoints=5000; % Calculating power in 5 sec. segments.
    SegmentN=ceil(length(LFP)/SegPoints);
    
    % calculate gamma and beta power
    parfor j=1:SegmentN-2
        
        local_LFP = LFP;
        
        segStart = (j-1)*SegPoints+1;
        segEnd = (j-1)*SegPoints+SegPoints;
        tempanaLFP = local_LFP(segStart:segEnd+dummypoints*2);
        
        % Low gamma.
        tempsignalLG = eegfilt(tempanaLFP, srate, 40, 65);
        tempLG = tempsignalLG(dummypoints+1:end-dummypoints); %%deleting the 2000 point in the beginning and the end
        signal_for_Lgamma = detrend(tempLG);
        PDdataLocal(j).logammaenergy = abs(hilbert(signal_for_Lgamma));
        PDdataLocal(j).logammaphase = angle(hilbert(signal_for_Lgamma));
        
        % High gamma.
        tempsignalHG = eegfilt(tempanaLFP,srate, 65, 100);
        tempHG = tempsignalHG(dummypoints+1:end-dummypoints); %%deleting the 2000 point in the beginning and the end
        signal_for_Hgamma = detrend(tempHG);
        PDdataLocal(j).higammaenergy = abs(hilbert(signal_for_Hgamma));
        PDdataLocal(j).higammaphase = angle(hilbert(signal_for_Hgamma));
        
        % Beta.
        tempsignalB = eegfilt(tempanaLFP,srate, 10,30);
        tempB = tempsignalB(dummypoints+1:end-dummypoints); %%%%deleting the 2000 point in the beginning and the end
        signal_for_beta = detrend(tempB);
        PDdataLocal(j).betaenergy = abs(hilbert(signal_for_beta));
        PDdataLocal(j).betaphase = angle(hilbert(signal_for_beta));
        
        % Delta.
        tempsignalD=eegfilt(tempanaLFP,srate, 1,4);
        tempD =tempsignalD(dummypoints+1:end-dummypoints); %%%%deleting the 2000 point in the beginning and the end
        signal_for_delta=detrend(tempD);
        PDdataLocal(j).deltaenergy = abs(hilbert(signal_for_delta));
        PDdataLocal(j).deltaphase= angle(hilbert(signal_for_delta));
        
        % Theta.
        tempsignalT = eegfilt(tempanaLFP,srate, 4,10);
        tempT = tempsignalT(dummypoints+1:end-dummypoints); %%%%deleting the 2000 point in the beginning and the end
        signal_for_theta = detrend(tempT);
        PDdataLocal(j).thetaenergy = abs(hilbert(signal_for_theta));
        PDdataLocal(j).thetaphase = angle(hilbert(signal_for_theta));

        % Low beta.
        tempsignalLB = eegfilt(tempanaLFP,srate, 10,20);
        tempLB = tempsignalLB(dummypoints+1:end-dummypoints); %%%%deleting the 2000 point in the beginning and the end
        signal_for_lobeta = detrend(tempLB);
        PDdataLocal(j).lobetaenergy = abs(hilbert(signal_for_lobeta));
        PDdataLocal(j).lobetaphase = angle(hilbert(signal_for_lobeta));
        
        % High beta.
        tempsignalHB = eegfilt(tempanaLFP,srate, 20,30);
        tempHB = tempsignalHB(dummypoints+1:end-dummypoints); %%%%deleting the 2000 point in the beginning and the end
        signal_for_hibeta = detrend(tempHB);
        PDdataLocal(j).hibetaenergy = abs(hilbert(signal_for_hibeta));
        PDdataLocal(j).hibetaphase = angle(hilbert(signal_for_hibeta));
        
        warning off all;
%         clear tempsignalHB tempsignalLB tempsignalB tempsignalLG tempsignalHG;
%         clear tempHB tempLB tempB tempLG tempHG;
%         clear signal_for_Lgamma signal_for_beta signal_for_Hgamma signal_for_lobeta signal_for_hibeta;
    end
    
    pdFields = fields(PDdataLocal);
    for n = 1:numel(pdFields)
        pdFieldName = pdFields{n};
        PDdata.(channel).(pdFieldName) = single(cat(2,PDdataLocal(:).(pdFieldName)));
    end
    
    % generates a spectrogram of the recording - Figure (24)
    fprintf('Generating a HT-Spectrogram for channel: %s\n',channel)
    PDdata.(channel).phase=[];
    PDdata.(channel).energy=[];
    
    [PDdata.(channel).phase,PDdata.(channel).energy] = basic_HT_improved_x11(LFP, srate, SegPoints, dummypoints, a, b, fileName);
    PDdata.(channel).phase = single(PDdata.(channel).phase);
    PDdata.(channel).energy = single(PDdata.(channel).energy);
    PDdata.(channel).LFP = single(LFP(dummypoints+1:end-dummypoints));
    
end


%     outputname=[DataAnalysisDirectory, fileName, '_HT', '.mat'];
    outputname = [fileName(1:end-4), '_HT.mat'];
    save (outputname, '-struct', 'PDdata', '-v6')
    
    cd (present_dir)
%     
%     clear phase; clear energy; clear ('betaphase', 'betaenergy', 'logammaenergy', 'logammaphase', 'higammaenergy', 'higammaphase','lobetaphase','lobetaenerg','hibetaphase','hibetaenergy');
% 
%     clear data;
% end