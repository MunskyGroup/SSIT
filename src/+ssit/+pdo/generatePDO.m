function [app,pdo] = generatePDO(app,paramsPDO,FSPoutputs,indsObserved,variablePDO)
arguments
    app
    paramsPDO = []
    FSPoutputs =[]
    indsObserved = []
    variablePDO = false
end
% This function generates the Probabilistic Distortion operator according
% to the user choices.

if isempty(FSPoutputs)
    if isempty(app.FspTabOutputs.solutions)
        app = runFsp(app);
    end
    FSPoutputs = app.FspTabOutputs.solutions;
end

%% Define Distortion Operator
nSpecies = ndims(FSPoutputs{1}.p.data);
if isempty(indsObserved)
    indsObserved = [1:nSpecies];
end
nObserved = length(indsObserved);

mxSize = zeros(1,nSpecies);
for it=length(FSPoutputs):-1:1
    for ispec = 1:nSpecies
        mxSize = max(mxSize,size(FSPoutputs{it}.p.data));
    end
end

conditionalPmfs = cell(1,nSpecies);

switch app.DistortionTypeDropDown.Value
    case 'None'
        for ispec = 1:nSpecies
            if mxSize(ispec)>1
                conditionalPmfs{ispec} = eye(mxSize(ispec));
            else
                conditionalPmfs{ispec} = 1;
            end
        end
    case 'Binomial'
        for ispec = 1:nSpecies
            if mxSize(ispec)>1
                lamb = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                for j = 1:mxSize(ispec)
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
                if mxSize(ispec)>1
                    lamb = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                    for j = 0:mxSize(ispec)-1
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
            if mxSize(ispec)>1
                lambFun = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                if ischar(lambFun)
                    lambFun = str2func(lambFun);
                end
                for j = 1:mxSize(ispec)
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
            if mxSize(ispec)>1
                lambFun = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                if ischar(lambFun)
                    lambFun = str2func(lambFun);
                end
                for j = 1:mxSize(ispec)
                    lamb = lambFun(j-1,paramsPDO);
                    conditionalPmfs{ispec}(1:j,j) = pdf('bino',0:j-1,j-1,lamb);
                end
            else
                conditionalPmfs{ispec} = 1;
            end
        end
    case 'Poisson'
        for ispec = 1:nSpecies
            if mxSize(ispec)>1
                lamb = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                for j = 1:mxSize(ispec)
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
                if mxSize(ispec)>1
                    lamb = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                    for j = 1:mxSize(ispec)
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
            if mxSize(ispec)>1
                lambFun = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                if ischar(lambFun)
                    lambFun = str2func(lambFun);
                end
                for j = 1:mxSize(ispec)
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
            if mxSize(ispec)>1
                lambFun = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                if ischar(lambFun)
                    lambFun = str2func(lambFun);
                end
                for j = 1:mxSize(ispec)
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
            if mxSize(ispec)>1
                nBins = app.FIMTabOutputs.PDOProperties.props.(['NumberBinsS',num2str(ispec)]);
                binType = app.FIMTabOutputs.PDOProperties.props.(['binTypeS',num2str(ispec)]);
                switch binType
                    case 'lin'
                        binEdges = unique(ceil(linspace(0,mxSize(ispec),nBins)));
                    case 'log'
                        binEdges = [0,unique(ceil(logspace(0,log10(mxSize(ispec)),nBins)))];
                end
                for j = 1:mxSize(ispec)
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
            if mxSize(ispec)>1
                conditionalPmfsBino=[];
                conditionalPmfsPoiss=[];
                dCdLamBino=[];
                dCdLamPoiss=[];
                alpha = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                for j = mxSize(ispec):-1:1
                    conditionalPmfsBino(1:j,j) = pdf('bino',0:j-1,j-1,alpha);
                end

                lamb = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                for j = mxSize(ispec):-1:1
                    numsFakeSpots = [0:ceil(lamb+5*sqrt(lamb))];
                    indsFakeSpots = j+numsFakeSpots;
                    probFakeSpots = pdf('poiss',numsFakeSpots,lamb);
                    conditionalPmfsPoiss(indsFakeSpots,j) = probFakeSpots;
                end
                conditionalPmfs{ispec} = conditionalPmfsPoiss*conditionalPmfsBino;

                if variablePDO
                    alpha = app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(ispec)]);
                    for j = mxSize(ispec)-1:-1:0
                        k=0:j;
                        dCdLamBino(k+1,j+1) = (k/alpha-(j-k)/(1-alpha)).*pdf('bino',k,j,alpha);
                    end
                    lamb = app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(ispec)]);
                    for j = mxSize(ispec):-1:1
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
            if mxSize(ispec)>1
                alpha = app.FIMTabOutputs.PDOProperties.props.PDOpars((ispec-1)*3+1);
                for j = mxSize(ispec):-1:1
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

    
    case 'Custom Function'
        for ispec = 1:nSpecies
            if mxSize(ispec)>1
                func = app.FIMTabOutputs.PDOProperties.props.(['PDO_S',num2str(ispec)]);
                if ischar(func)
                    func = str2func(func);
                end
                maxObs = app.FIMTabOutputs.PDOProperties.props.(['MaxObservationS',num2str(ispec)]);
                for j = 1:mxSize(ispec)
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