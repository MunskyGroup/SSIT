function output = Simulate(kon, koff, w, kex, kr, D, gam, NCells, makePlot)
    a = 53*ones(NCells,1);
    b = 39*ones(NCells,1);
    u = 2*pi*rand(NCells,1);
    posnTS=[a.*cos(u),b.*sin(u)].*rand(NCells,2); %position of TS for sim data

    outputFolder = 'SimOutput';
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Generate a unique timestamp string with milliseconds
    timestamp = datetime("now", 'yyyymmdd_HHMMSSFFF');
    
    % Generate a random number to reduce the chance of collision
    randomNumber = randi([1000, 9999]);
    
    % Construct the unique filename
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
