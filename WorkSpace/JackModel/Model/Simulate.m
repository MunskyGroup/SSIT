function output = Simulate(kon, koff, w, kex, kr, D, gam, NCells, makePlot)
    a = 53*ones(NCells,1);
    b = 39*ones(NCells,1);
    u = 2*pi*rand(NCells,1);
    posnTS=[a.*cos(u),b.*sin(u)].*rand(NCells,2); %position of TS for sim data

    outputFolder = 'SimOutput';
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Generate timestamp with milliseconds
    timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmssSSS'));
    
    % Generate random 4-digit number
    randomNumber = randi([1000, 9999]);
    
    % Create unique filename
    fileName = sprintf('Data_%s_%d.csv', timestamp, randomNumber);
    
    % Full file path
    fullFilePath = fullfile(outputFolder, fileName);


    stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,NCells,fullFilePath)
    
    [intensTS,binCounts,abEllipse,areaEllipse,ecentricity,...
        binTS,binCountsDim,binCountsBright,hisCountsBright,...
        hisCountsDim,XcRot,iTS,XCVcRot,XYfit] = ...
        processSpotPositionData(fullFilePath,[],a,b,0,11,0);
    
    output.intensTS = intensTS;
    output.binCounts = binCounts;
    output.abEllipse = abEllipse;
    output.areaEllipse = areaEllipse;
    output.ecentricity = ecentricity;
    output.binTS = binTS;
    output.binCountsDim = binCountsDim;
    output.binCountsBright = binCountsBright;
    output.hisCountsBright = hisCountsBright;
    output.hisCountsDim = hisCountsDim;
    output.XcRot = XcRot;
    output.iTS = iTS;
    output.XCVcRot = XCVcRot;
    output.XYfit = XYfit;
    output.filePath = fullFilePath;
end
