function output = Simulate(kon, koff, w, kex, kr, D, gam, NCells, makePlot)
    a = 53*ones(NCells,1);
    b = 39*ones(NCells,1);
    u = 2*pi*rand(NCells,1);
    posnTS=[a.*cos(u),b.*sin(u)].*rand(NCells,2); %position of TS for sim data

    outputFolder = 'SimOutput';
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % List all existing files in the folder that match "Data*.csv"
    filePattern = fullfile(outputFolder, 'Data*.csv');
    existingFiles = dir(filePattern);
    
    % Extract numbers from filenames like "Data12.csv"
    existingNums = zeros(length(existingFiles), 1);
    for i = 1:length(existingFiles)
        name = existingFiles(i).name;
        tokens = regexp(name, 'Data(\d+)\.csv', 'tokens');
        if ~isempty(tokens)
            existingNums(i) = str2double(tokens{1}{1});
        end
    end
    
    % Determine the next available number
    if isempty(existingNums)
        nextNum = 1;
    else
        nextNum = max(existingNums) + 1;
    end
    
    % Create filename
    fileName = sprintf('Data%d.csv', nextNum);
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
