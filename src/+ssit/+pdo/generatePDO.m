function [app,pdo] = generatePDO(app,paramsPDO,FSPoutputs,indsObserved,variablePDO,maxSize)
arguments
    app
    paramsPDO = []
    FSPoutputs =[]
    indsObserved = []
    variablePDO = false
    maxSize=[]
end
% This function generates the Probabilistic Distortion operator according
% to the user choices.

if isempty(maxSize)
    if isempty(FSPoutputs)
        if isempty(app.FspTabOutputs.solutions)
            app = runFsp(app);
        end
        FSPoutputs = app.FspTabOutputs.solutions;
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
                lamb = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
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
                    lamb = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
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
if variablePDO
    app.FIMTabOutputs.distortionOperator = ssit.pdo.TensorProductDistortionOperator(conditionalPmfs(indsObserved),dCdLam(indsObserved,:));
else
    app.FIMTabOutputs.distortionOperator = ssit.pdo.TensorProductDistortionOperator(conditionalPmfs(indsObserved),[]);
end

if nargout>1
    pdo = app.FIMTabOutputs.distortionOperator;
end