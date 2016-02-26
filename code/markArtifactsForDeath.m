% some code borrowed from https://gist.github.com/cholland29/3107790
function [artifacts] = markArtifactsForDeath(X, samplingFreq, segmentTime)
    artifacts = struct('indices',{},...
                    'm1',{},...
                    'str',{}...
                    );
    % following function borrowed from https://gist.github.com/cholland29/3107790
    function OnClickAxes(hax, evt)
        point1 = get(hax,'CurrentPoint'); % corner where rectangle starts ( initial mouse down point )
        rbbox
        point2 = get(hax,'CurrentPoint'); % corner where rectangle stops ( when user lets go of mouse )

        % Now lets iterate through all lines in the axes and extract the data that lies within the selected region
        allLines = findall(hax,'type','line');
        
        dataInd = getDataInRect( point1(1,1:2), point2(1,1:2), allLines(1) ); % not interested in z-coord
        
        tempCell = {allLines.YData};
        shapeMat = cell2mat(tempCell(:));
        artifacts(end+1).indices = dataInd;
        artifacts(end).m1 = shapeMat(2,dataInd);
        artifacts(end).str = shapeMat(1,dataInd);
    end
    
    function dataInd = getDataInRect(point1, point2, line)
        start = min(point1(1),point2(1));

        fin = max(point1(1), point2(1));

        dataInd = find((line.XData >= start) & (line.XData <= fin));
        fprintf('Range from %.0f to %.0f added\n',start,fin);
    end


    t = samplingFreq*(1:1:length(X));

    if isempty(segmentTime)
        segmentTime = 30; % seconds
    end

    segmentIndices = segmentTime*samplingFreq; 
    h = figure;
    hax = axes;
    plt = plot(t,X); xlabel('Time [s]'), ylabel('Voltage [mV]')
    set(hax, 'ButtonDownFcn',@OnClickAxes);
    pause
end



