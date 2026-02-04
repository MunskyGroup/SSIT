function [app,pdo] = generatePDO(app,paramsPDO,FSPoutputs,indsObserved,variablePDO,maxSize,showPlot,opts)
%% SSIT.generatePDO - This function generates the Probabilistic 
%% Distortion Operator (PDO) according to user choice.
% 
% Inputs:
%   * app
%   * paramsPDO - ()
%   * FSPoutputs - ()
%   * indsObserved - ()
%   * variablePDO - (logical), default: false
%   * maxSize - ()
%   * showPlot - (logical), default: true
%
%   Optional plotting arguments:
%   * opts.Title - (string) 
%   * opts.FontSize - (double), default: 18
%   * opts.XLabel - (string), default: "True counts"
%   * opts.YLabel - (string), default: "Observed counts"
%   * opts.XLim - (double), default: []
%   * opts.YLim - (double), default: []
%   * opts.Levels - (double), default: 30
%   * opts.CLim double - default: [-25 0] 
% Output: 
%
% Example: 
arguments
    app
    paramsPDO = []
    FSPoutputs = []
    indsObserved = []
    variablePDO = false
    maxSize = []
    showPlot (1,1) logical = true

    % Plotting opts (used only when showPlot==true)
    opts.Title (1,1) string = ""
    opts.FontSize (1,1) double {mustBePositive} = 18
    opts.XLabel (1,1) string = "True counts"
    opts.YLabel (1,1) string = "Observed counts"
    opts.XLim double = []
    opts.YLim double = []
    opts.Levels (1,1) double {mustBeInteger, mustBePositive} = 30
    opts.CLim double = [-25 0]   % log10(P) clamp range for display
end


if isempty(maxSize)
    if isempty(FSPoutputs)
        if isempty(app.SSITModel.Solutions)
            app = runFsp(app);
        end
        FSPoutputs = app.SSITModel.Solutions.fsp;
    end
    nSpecies = ndims(FSPoutputs{1}.p.data);
    maxSize = zeros(1,nSpecies);
    for it=length(FSPoutputs):-1:1
        for ispec = 1:nSpecies
            maxSize = max(maxSize,size(FSPoutputs{it}.p.data));
        end
    end
else
    nSpecies = length(maxSize);
end

%% Define Distortion Operator
if isempty(indsObserved)
    indsObserved = [1:nSpecies];
end
% nObserved = length(indsObserved);

conditionalPmfs = cell(1,nSpecies);

switch app.DistortionTypeDropDown.Value
    case 'None'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                conditionalPmfs{ispec} = eye(maxSize(ispec));
            else
                conditionalPmfs{ispec} = 1;
            end
        end
    case 'Binomial'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                % --- Prefer calibrated parameters if provided ---
                if ~isempty(paramsPDO)
                    lamb = paramsPDO(ispec);   % lamb = capture probability p_cap for species ispec
                elseif isfield(app.FIMTabOutputs.PDOProperties.props,'PDOpars') && ...
                       numel(app.FIMTabOutputs.PDOProperties.props.PDOpars) >= ispec
                    lamb = app.FIMTabOutputs.PDOProperties.props.PDOpars(ispec);
                else
                    lamb = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                end
                lamb = min(max(lamb,0),1); % clamp just in case
                % ------------------------------------------------------                
                for j = 1:maxSize(ispec)
                    conditionalPmfs{ispec}(1:j,j) = pdf('bino',0:j-1,j-1,lamb);
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end

        if variablePDO
            dCdLam = cell(nSpecies,nSpecies);
            for ispec = 1:nSpecies
                for jspec = 1:nSpecies
                    dCdLam{jspec,ispec} = [];
                end
                if maxSize(ispec)>1                    
                    % --- Prefer calibrated parameters if provided ---
                    if ~isempty(paramsPDO)
                        lamb = paramsPDO(ispec);   % lamb = capture probability p_cap for species ispec
                    elseif isfield(app.FIMTabOutputs.PDOProperties.props,'PDOpars') && ...
                        numel(app.FIMTabOutputs.PDOProperties.props.PDOpars) >= ispec
                        lamb = app.FIMTabOutputs.PDOProperties.props.PDOpars(ispec);
                    else
                        lamb = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                    end
                lamb = min(max(lamb,0),1); % clamp just in case
                % ------------------------------------------------------ 
                    for j = 0:maxSize(ispec)-1
                        pdf('bino',0:j-1,j-1,lamb);
                        k=0:j;
                        dCdLam{ispec,ispec}(k+1,j+1) = (k/lamb-(j-k)/(1-lamb)).*pdf('bino',k,j,lamb);
                    end
                else
                    dCdLam{ispec,ispec} = 0;
                end
            end
        end

    case 'Binomial - State Dependent'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                lambFun = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                if ischar(lambFun)
                    lambFun = str2func(lambFun);
                end
                for j = 1:maxSize(ispec)
                    lamb = lambFun(j-1);
                    conditionalPmfs{ispec}(1:j,j) = pdf('bino',0:j-1,j-1,lamb);
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end
    case 'Binomial - Parametrized'
        if isempty(paramsPDO)
            if isnumeric(app.FIMTabOutputs.PDOProperties.props.ParameterGuess)
                paramsPDO = app.FIMTabOutputs.PDOProperties.props.ParameterGuess;
            else
                paramsPDO = eval(app.FIMTabOutputs.PDOProperties.props.ParameterGuess);
            end
        end
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                lambFun = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                if ischar(lambFun)
                    lambFun = str2func(lambFun);
                end
                for j = 1:maxSize(ispec)
                    if isnumeric(lambFun)
                        lamb = lambFun;
                    else
                        lamb = lambFun(j-1,paramsPDO);
                    end
                    conditionalPmfs{ispec}(1:j,j) = pdf('bino',0:j-1,j-1,lamb);
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end
    case 'Poisson'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                lamb = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                for j = 1:maxSize(ispec)
                    numsFakeSpots = [0:ceil(lamb+5*sqrt(lamb))];
                    indsFakeSpots = j+numsFakeSpots;
                    probFakeSpots = pdf('poiss',numsFakeSpots,lamb);
                    conditionalPmfs{ispec}(indsFakeSpots,j) = probFakeSpots;
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end
        if variablePDO
            dCdLam = cell(nSpecies,nSpecies);
            for ispec = 1:nSpecies
                for jspec = 1:nSpecies
                    dCdLam{jspec,ispec} = [];
                end
                if maxSize(ispec)>1
                    lamb = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                    for j = 1:maxSize(ispec)
                        numsFakeSpots = [0:ceil(lamb+5*sqrt(lamb))];
                        indsFakeSpots = j+numsFakeSpots;
                        probFakeSpots = pdf('poiss',numsFakeSpots,lamb);
                        dCdLam{ispec,ispec}(indsFakeSpots,j) = (numsFakeSpots/lamb-1).*probFakeSpots;
                    end
                else
                    dCdLam{ispec,ispec} = 0;
                end
            end
        end
    case 'Poisson - State Dependent'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                lambFun = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                if ischar(lambFun)
                    lambFun = str2func(lambFun);
                end
                for j = 1:maxSize(ispec)
                    lamb = lambFun(j-1);
                    numsFakeSpots = [0:ceil(lamb+5*sqrt(lamb))];
                    indsFakeSpots = j+numsFakeSpots;
                    probFakeSpots = pdf('poiss',numsFakeSpots,lamb);
                    conditionalPmfs{ispec}(indsFakeSpots,j) = probFakeSpots;
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end
    case 'Poisson - Parametrized'
        if isempty(paramsPDO)
            if isnumeric(app.FIMTabOutputs.PDOProperties.props.ParameterGuess)
                paramsPDO = app.FIMTabOutputs.PDOProperties.props.ParameterGuess;
            else
                paramsPDO = eval(app.FIMTabOutputs.PDOProperties.props.ParameterGuess);
            end
        end
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                lambFun = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                if ischar(lambFun)
                    lambFun = str2func(lambFun);
                end
                for j = 1:maxSize(ispec)
                    lamb = lambFun(j-1,paramsPDO);
                    numsFakeSpots = [0:ceil(lamb+5*sqrt(lamb))];
                    indsFakeSpots = j+numsFakeSpots;
                    probFakeSpots = pdf('poiss',numsFakeSpots,lamb);
                    conditionalPmfs{ispec}(indsFakeSpots,j) = probFakeSpots;
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end
    case 'Binning'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                nBins = app.FIMTabOutputs.PDOProperties.props.(['NumberBinsS',num2str(ispec)]);
                binType = app.FIMTabOutputs.PDOProperties.props.(['binTypeS',num2str(ispec)]);
                switch binType
                    case 'lin'
                        binEdges = unique(ceil(linspace(0,maxSize(ispec),nBins)));
                    case 'log'
                        binEdges = [0,unique(ceil(logspace(0,log10(maxSize(ispec)),nBins)))];
                end
                for j = 1:maxSize(ispec)
                    k = find(j>=binEdges,1,"last");
                    conditionalPmfs{ispec}(k,j) = 1;
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end

    case 'BinoPoiss'
        if variablePDO
            dCdLam = cell(nSpecies,2*nSpecies);
            for ispec = 1:nSpecies
                for jspec = 1:nSpecies
                    dCdLam{ispec,jspec} = [];
                end
            end
        end
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                conditionalPmfsBino=[];
                conditionalPmfsPoiss=[];
                dCdLamBino=[];
                dCdLamPoiss=[];
                alpha = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                for j = maxSize(ispec):-1:1
                    conditionalPmfsBino(1:j,j) = pdf('bino',0:j-1,j-1,alpha);
                end

                lamb = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                for j = maxSize(ispec):-1:1
                    numsFakeSpots = [0:ceil(lamb+5*sqrt(lamb))];
                    indsFakeSpots = j+numsFakeSpots;
                    probFakeSpots = pdf('poiss',numsFakeSpots,lamb);
                    conditionalPmfsPoiss(indsFakeSpots,j) = probFakeSpots;
                end
                conditionalPmfs{ispec} = conditionalPmfsPoiss*conditionalPmfsBino;

                if variablePDO
                    alpha = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                    for j = maxSize(ispec)-1:-1:0
                        k=0:j;
                        dCdLamBino(k+1,j+1) = (k/alpha-(j-k)/(1-alpha)).*pdf('bino',k,j,alpha);
                    end
                    lamb = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                    for j = maxSize(ispec):-1:1
                        numsFakeSpots = [0:ceil(lamb+5*sqrt(lamb))];
                        indsFakeSpots = j+numsFakeSpots;
                        probFakeSpots = pdf('poiss',numsFakeSpots,lamb);
                        dCdLamPoiss(indsFakeSpots,j) = (numsFakeSpots/lamb-1).*probFakeSpots;
                    end
                    dCdLam{ispec,(ispec-1)*2+1} =  conditionalPmfsPoiss*dCdLamBino;
                    dCdLam{ispec,ispec*2} =  dCdLamPoiss*conditionalPmfsBino;
                end
            else
                conditionalPmfs{ispec} = 1;
                if variablePDO
                    dCdLam{ispec,ispec*2-1} = 0;
                    dCdLam{ispec,ispec*2} = 0;
                end
            end
        end


    case 'BinoAffinePoiss'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                alpha = app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+1);
                for j = maxSize(ispec):-1:1
                    lamb = app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+2)+...
                        app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+3)*(j-1);
                    Np = ceil(lamb+10*sqrt(lamb))+50;
                    P2 = pdf('poiss',[0:Np],lamb);
                    P1 = binopdf([0:j-1],j-1,alpha);
                    conditionalPmfs{ispec}(1:j+Np,j) = conv(P2,P1);
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end
    case 'AffinePoissLoss'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                alpha = app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+3);
                for j = maxSize(ispec):-1:1
                    lamb = app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+1)+...
                        app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+2)*(j-1);
                    Np = ceil(lamb+10*sqrt(lamb))+50;
                    P2 = pdf('poiss',[0:Np],lamb);
                    conditionalPmfs{ispec}(1:Np+1,j) = P2*alpha;
                    conditionalPmfs{ispec}(1,j) = conditionalPmfs{ispec}(1,j) + 1-alpha;
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end

    case 'DiscretizedNormal'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                lambda = app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*4+[1:4]);
                for j = maxSize(ispec):-1:1
                    mu = lambda(1)+lambda(3)*(j-1);
                    sig = lambda(2)+lambda(4)*(j-1)+1e-5;
                    if isfield(app.FIMTabOutputs.PDOProperties.props,'pdoOutputRange')
                        Np = app.FIMTabOutputs.PDOProperties.props.pdoOutputRange(ispec);
                    else
                        Np = ceil(mu+10*sqrt(sig))+50;
                    end
                    P2 = cdf("Normal",[0:Np]+.5,mu,sqrt(sig))-cdf("Normal",[0:Np]-.5,mu,sqrt(sig));
%                     P2 = pdf('normal',[0:Np],mu,sqrt(sig2));
                    conditionalPmfs{ispec}(1:Np+1,j) = P2;
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end

    case 'GaussSpurrious'
        NFail = 10;
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                lambda = app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*7+[1:7]);
                for j = maxSize(ispec):-1:1
                    mu = lambda(1)+lambda(3)*(j-1);
                    sig = sqrt(lambda(2)+lambda(4)*(j-1));
                    muFail = lambda(6);
                    sigFail = sqrt(lambda(7));
                    if isfield(app.FIMTabOutputs.PDOProperties.props,'pdoOutputRange')
                        Np = app.FIMTabOutputs.PDOProperties.props.pdoOutputRange(ispec);
                    else
                        Np = ceil(mu+10*sqrt(sig))+50;
                    end
                    P = zeros(Np+1,1);
                    edges = [0:max(9,Np)];
                    if j<=NFail
                        P(1) = (1-lambda(5))*(cdf('norm',.5,mu,sig))+...
                            lambda(5)*(cdf('norm',.5,muFail,sigFail));
                        P(2:9) = (1-lambda(5))*(cdf('norm',edges(2:9)+.5,mu,sig)-cdf('norm',edges(2:9)-.5,mu,sig))+...
                            lambda(5)*(cdf('norm',edges(2:9)+.5,muFail,sigFail)-cdf('norm',edges(2:9)-.5,muFail,sigFail));
                        P(10:end) = (1-lambda(5))*pdf('norm',edges(10:end),mu,sig)+...
                            lambda(5)*pdf('norm',edges(10:end),muFail,sigFail);
                    else
                        P(1) = cdf('norm',.5,mu,sig);
                        P(2:9) = cdf('norm',edges(2:9)+.5,mu,sig)-cdf('norm',edges(2:9)-.5,mu,sig);
                        P(10:end) = pdf('norm',edges(10:end),mu,sig);
                    end
                    conditionalPmfs{ispec}(1:Np+1,j) = P(1:Np+1);
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end

    case 'AffinePoiss'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                for j = maxSize(ispec):-1:1
                    lamb = max(app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+1),...
                        app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+2)+...
                        app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+3)*(j-1));
                    if isfield(app.FIMTabOutputs.PDOProperties.props,'pdoOutputRange')
                        Np = app.FIMTabOutputs.PDOProperties.props.pdoOutputRange(ispec);
                    else
                        Np = ceil(lamb+10*sqrt(lamb))+50;
                    end
                    P2 = pdf('poiss',[0:Np],lamb);
                    conditionalPmfs{ispec}(1:Np+1,j) = P2;
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end


    case 'AffinePoissBounded'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                for j = maxSize(ispec):-1:1
                    lamb = max(1,app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*2+1)+...
                        app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*2+2)*(j-1));
                    Np = ceil(lamb+10*sqrt(lamb))+50;
                    P2 = pdf('poiss',[0:Np],lamb);
                    conditionalPmfs{ispec}(1:Np+1,j) = P2;
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end


    case 'Custom Function'
        for ispec = 1:nSpecies
            if maxSize(ispec)>1
                func = app.FIMTabOutputs.PDOProperties.props.(['PDO_S',num2str(ispec)]);
                if ischar(func)
                    func = str2func(func);
                end
                maxObs = app.FIMTabOutputs.PDOProperties.props.(['MaxObservationS',num2str(ispec)]);
                for j = 1:maxSize(ispec)
                    for i = 0:maxObs
                        try
                            conditionalPmfs{ispec}(i+1,j) = func(i,j-1);
                        catch
                            conditionalPmfs{ispec}(i+1,j) = 0;
                        end
                    end
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end
end
% --------------------------------------------------------------------
% Optional visualization of per-species conditional PMFs
%   conditionalPmfs{ispec} is (nObserved x nTrue)
%   x-axis = true counts (columns), y-axis = observed counts (rows)
% --------------------------------------------------------------------
if showPlot
    if isempty(indsObserved)
        indsToPlot = 1:nSpecies;
    else
        indsToPlot = indsObserved;
    end

    % Only plot species with non-trivial matrices
    indsToPlot = indsToPlot(:)';
    indsToPlot = indsToPlot(indsToPlot >= 1 & indsToPlot <= nSpecies);

    % Make one figure with a tiled layout
    fg = figure;
    tl = tiledlayout(fg, 'flow');
    if strlength(opts.Title) > 0
        title(tl, opts.Title, 'FontSize', opts.FontSize);
    end

    for ii = 1:numel(indsToPlot)
        ispec = indsToPlot(ii);
        P = conditionalPmfs{ispec};

        % Skip scalar/singleton PDOs
        if isempty(P) || (isscalar(P) && numel(P)==1)
            nexttile;
            axis off;
            text(0.1, 0.5, sprintf('S%d: (trivial PDO)', ispec), 'FontSize', opts.FontSize);
            continue
        end

        % Log10 display with safe floor to avoid -Inf
        Z = log10(max(P, realmin));
        Z = max(opts.CLim(1), min(opts.CLim(2), Z));  % clamp for nicer contrast

        nexttile;
        contourf(0:size(P,2)-1, 0:size(P,1)-1, Z, opts.Levels, 'LineColor','none');
        colorbar;
        xlabel(opts.XLabel);
        ylabel(opts.YLabel);

        if ~isempty(opts.XLim); xlim(opts.XLim); end
        if ~isempty(opts.YLim); ylim(opts.YLim); end

        set(gca, 'FontSize', max(10, round(0.8*opts.FontSize)));

        % Per-panel title
        if strlength(opts.Title) == 0
            title(sprintf('PDO (S%d)', ispec), 'FontSize', opts.FontSize);
        else
            title(sprintf('S%d', ispec), 'FontSize', opts.FontSize);
        end
    end
end
if variablePDO
    app.FIMTabOutputs.distortionOperator = ssit.pdo.TensorProductDistortionOperator(conditionalPmfs(indsObserved),dCdLam(indsObserved,:));
else
    app.FIMTabOutputs.distortionOperator = ssit.pdo.TensorProductDistortionOperator(conditionalPmfs(indsObserved),[]);
end

if nargout>1
    pdo = app.FIMTabOutputs.distortionOperator;
end