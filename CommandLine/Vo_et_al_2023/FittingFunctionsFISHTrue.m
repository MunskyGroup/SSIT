function results = FittingFunctionsFISHTrue(iStep,fileName,vars)
arguments
    iStep
    fileName
    vars = [];
end

try
    copyfile([fileName,'.mat'],[fileName,'_BAK.mat']);
catch
    disp('Previous file not available -- will need to start from scratch.')
end

results = [];

defaults = {'priorScale',1;...
    'timeSet',[0,18,300];...
    'nMH',3000;...
    'nFIMsamples',4;...
    'modelVarsToFit',[1:5];...
    'iter',400;...
    'display','final';...
    'mleScatterPlots',true;...
    'mhScaling',0.6;...
    'pdoTimes',[0,300];...
    'FIMnCells',[0,0,0];...
    'dataFile','Huy_intensity_data_correct/NucAndSpotClassification_dTS_2.csv';...
    'modelChoice','2stateBurst';...
    'muLog10Prior',[-4,-4,log10(0.2),log10(7.12),log10(0.012)]';...
    'sigLog10Prior',[1,1,.5,.5,.5]'*vars.priorScale;...
    'fitCases',[1:5];...
    'makePlots',false;...
    'otherParGuesses',[];...
    'initialFspBounds',[];...
    'covWithPrior',false};

for i=1:size(defaults,1)
    if ~isfield(vars,defaults{i,1})
        vars.(defaults{i,1}) = defaults{i,2};
    end
end
binSize = 20;
chainResults=cell(1,5);
fitOptions = optimset('Display',vars.display,'MaxIter',vars.iter);
muLog10Prior=vars.muLog10Prior;
sigLog10Prior=vars.sigLog10Prior;

%%
switch iStep
    case 1 % Fit the PDO parameters to real and distorted data.
        try
            load(fileName,'pdoPars')
        catch
            warning('Previous PDO fits not found -- restarting.')
            pdoPars.lPoissTrue2MCPspots = [1,0.5,1.0];
            pdoPars.lPoissTrue2FISHintens = [1,1000,1.0];
            pdoPars.lPoissTrue2MCPintens = [1,1000,1.0];
            pdoPars.lNormTrue2MCPspots = [1,1,1,1];
            pdoPars.lNormTrue2FISHintens = [1000,1e6,1,1];
            pdoPars.lNormTrue2MCPintens = [1000,1e6,1,1];
        end

        if ~isfield(pdoPars,'lNormSpurriousMCPspots')
            pdoPars.lNormSpurriousMCPspots = [1,1,1,1,0.01,10,100];
        end        
        if ~isfield(pdoPars,'lNormSpurriousMCPintens')
            pdoPars.lNormSpurriousMCPintens = [1000,1e6,1,1,0.01,1000,1e6];
        end        
        if ~isfield(pdoPars,'lNormSpurriousFISHintens')
            pdoPars.lNormSpurriousFISHintens = [1000,1e6,1,1,0.01,1000,1e6];
        end        
        if ~isfield(pdoPars,'lNormAddMCPspots')
            pdoPars.lNormAddMCPspots = [0.8,25,1000];
        end        

        DATA = importdata(vars.dataFile);
        % 0 = MCP-GFP
        % 1 = smiFISH
        K = find(strcmp(DATA.colheaders,'number_spots_type_0')); % Column with 'MCP-GFP' data
        K2 = find(strcmp(DATA.colheaders,'nucIntens1')); % Column with 'MCP' Intensity data

        J2 = find(strcmp(DATA.colheaders,'nucIntens3')); % Column with 'smiFISH' Intensity data

        I = find(strcmp(DATA.colheaders,'number_spots_type_1'));  % Column with 'true' data

        jT = find(strcmp(DATA.colheaders,'time'));

        jCells = DATA.data(:,jT)==vars.pdoTimes(1);
        for i=2:length(vars.pdoTimes)
            jCells(DATA.data(:,jT)==vars.pdoTimes(i))=true;
        end

        DistortedFISHintens = DATA.data(jCells,J2);
        DistortedMCPspots = DATA.data(jCells,K);
        DistortedMCPintens = DATA.data(jCells,K2);
        True = DATA.data(jCells,I);

        limsMcp = [0,max(DATA.data(:,K))];
        limsIfish = [0,max(DATA.data(:,J2))];
        limsImcp = [0,max(DATA.data(:,K2))];

        if vars.doFit==1
            % i = FISH,MCP
            % j = spots,intens
            % k = poiss,norm
            for iFit = 1:10
                OBJ = @(l)-findError(l,True,DistortedFISHintens,limsIfish);
                pdoPars.lPoissTrue2FISHintens = fminsearch(OBJ,pdoPars.lPoissTrue2FISHintens,fitOptions);
                OBJ = @(l)-findErrorNorm(l,True,DistortedFISHintens,limsIfish);
                pdoPars.lNormTrue2FISHintens = fminsearch(OBJ,pdoPars.lNormTrue2FISHintens,fitOptions);
                OBJ = @(l)-findErrorNormSpurrious(l,True,DistortedFISHintens,limsIfish);
                pdoPars.lNormSpurriousFISHintens = fminsearch(OBJ,pdoPars.lNormSpurriousFISHintens,fitOptions);
                
                OBJ = @(l)-findError(l,True,DistortedMCPspots,limsMcp);
                pdoPars.lPoissTrue2MCPspots = fminsearch(OBJ,pdoPars.lPoissTrue2MCPspots,fitOptions);
                OBJ = @(l)-findErrorNorm(l,True,DistortedMCPspots,limsMcp);
                pdoPars.lNormTrue2MCPspots = fminsearch(OBJ,pdoPars.lNormTrue2MCPspots,fitOptions);
                OBJ = @(l)-findErrorNormSpurrious(l,True,DistortedMCPspots,limsMcp);
                pdoPars.lNormSpurriousMCPspots = fminsearch(OBJ,pdoPars.lNormSpurriousMCPspots,fitOptions);
                OBJ = @(l)-findErrorNormAdd(l,True,DistortedMCPspots,limsMcp);
                pdoPars.lNormAddMCPspots = fminsearch(OBJ,pdoPars.lNormAddMCPspots,fitOptions);

                OBJ = @(l)-findError(l,True,DistortedMCPintens,limsImcp);
                pdoPars.lPoissTrue2MCPintens = fminsearch(OBJ,pdoPars.lPoissTrue2MCPintens,fitOptions);
                OBJ = @(l)-findErrorNorm(l,True,DistortedMCPintens,limsImcp);
                pdoPars.lNormTrue2MCPintens = fminsearch(OBJ,pdoPars.lNormTrue2MCPintens,fitOptions);
                OBJ = @(l)-findErrorNormSpurrious(l,True,DistortedMCPintens,limsIfish);
                pdoPars.lNormSpurriousMCPintens = fminsearch(OBJ,pdoPars.lNormSpurriousMCPintens,fitOptions);
            end
            try
                save(fileName,'pdoPars','-append')
            catch
                save(fileName,'pdoPars');
            end
            pdoPars
        else

            [logLFit(1),cSet{1}] = findError(pdoPars.lPoissTrue2MCPspots,True,DistortedMCPspots,2*limsMcp);
            [logLFit(2),cSet{2}] = findErrorNorm(pdoPars.lNormTrue2MCPspots,True,DistortedMCPspots,2*limsMcp);
            [logLFit(3),cSet{3}] = findErrorNormSpurrious(pdoPars.lNormSpurriousMCPspots,True,DistortedMCPspots,2*limsMcp);
            [logLFit(10),cSet{10}] = findErrorNormAdd(pdoPars.lNormAddMCPspots,True,DistortedMCPspots,2*limsMcp);
           
            [logLFit(4),cSet{4}] = findError(pdoPars.lPoissTrue2FISHintens,True,DistortedFISHintens,2*limsIfish);
            [logLFit(5),cSet{5}] = findErrorNorm(pdoPars.lNormTrue2FISHintens,True,DistortedFISHintens,2*limsIfish);
            [logLFit(6),cSet{6}] = findErrorNormSpurrious(pdoPars.lNormSpurriousFISHintens,True,DistortedFISHintens,2*limsIfish);

            [logLFit(7),cSet{7}] = findError(pdoPars.lPoissTrue2MCPintens,True,DistortedMCPintens,2*limsImcp);
            [logLFit(8),cSet{8}] = findErrorNorm(pdoPars.lNormTrue2MCPintens,True,DistortedMCPintens,2*limsImcp);
            [logLFit(9),cSet{9}] = findErrorNormSpurrious(pdoPars.lNormSpurriousMCPintens,True,DistortedMCPintens,2*limsImcp);

            results.loglFit = logLFit;

            close all
            logLPred = zeros(10,3);
            for iTime = 1:3
                switch iTime
                    case 1
                        jCells = DATA.data(:,jT)==0;
                    case 2
                        jCells = DATA.data(:,jT)==18;
                    case 3
                        jCells = DATA.data(:,jT)==300;
                end

                DistortedFISHintens = DATA.data(jCells,J2);
                DistortedMCPspots = DATA.data(jCells,K);
                DistortedMCPintens = DATA.data(jCells,K2);
                True = DATA.data(jCells,I);

                [logLPred(1,iTime)] = findError(pdoPars.lPoissTrue2MCPspots,True,DistortedMCPspots,2*limsMcp);
                [logLPred(2,iTime)] = findErrorNorm(pdoPars.lNormTrue2MCPspots,True,DistortedMCPspots,2*limsMcp);
                [logLPred(3,iTime)] = findErrorNormSpurrious(pdoPars.lNormSpurriousMCPspots,True,DistortedMCPspots,2*limsMcp);
                [logLPred(10,iTime)] = findErrorNormAdd(pdoPars.lNormAddMCPspots,True,DistortedMCPspots,2*limsMcp);

                [logLPred(4,iTime)] = findError(pdoPars.lPoissTrue2FISHintens,True,DistortedFISHintens,2*limsIfish);
                [logLPred(5,iTime)] = findErrorNorm(pdoPars.lNormTrue2FISHintens,True,DistortedFISHintens,2*limsIfish);
                [logLPred(6,iTime)] = findErrorNormSpurrious(pdoPars.lNormSpurriousFISHintens,True,DistortedFISHintens,2*limsIfish);

                [logLPred(7,iTime)] = findError(pdoPars.lPoissTrue2MCPintens,True,DistortedMCPintens,2*limsImcp);
                [logLPred(8,iTime)] = findErrorNorm(pdoPars.lNormTrue2MCPintens,True,DistortedMCPintens,2*limsImcp);
                [logLPred(9,iTime)] = findErrorNormSpurrious(pdoPars.lNormSpurriousMCPintens,True,DistortedMCPintens,2*limsImcp);

                results.loglPred = logLPred;

                if vars.makePlots
                    DistSet = {DistortedMCPspots,DistortedMCPspots,DistortedMCPspots,...
                        DistortedFISHintens,DistortedFISHintens,DistortedFISHintens,...
                        DistortedMCPintens,DistortedMCPintens,DistortedMCPintens,DistortedMCPspots};
                    labSet = {'MCP-GFP','MCP-GFP','MCP-GFP',...
                        'I_{FISH}','I_{FISH}','I_{FISH}',...
                        'I_{MCP-GFP}','I_{MCP-GFP}','I_{MCP-GFP}','MCP-GFP'};
                    titles = {'MCP-GFP Spots, Poisson','MCP-GFP Spots, Gaussian','MCP-GFP Spots, Gauss+Spurrious',...
                        'FISH Intensity, Poisson','FISH Intensity, Gaussian','FISH Intensity, Gaussian+Spurrious',...
                        'MCP-GFP Intensity, Poisson','MCP-GFP Intensity, Gaussian','MCP-GFP Intensity, Gaussian+Spurrious','MCP-GFP Spots, NormAdd'};

                    col=[0.2,0.2,0.6;...
                        0.6,0.2,0.2;...
                        0.2,0.6,0.6;...
                        0.4,0.4,0.6;...
                        0.6,0.4,0.4;...
                        0.4,0.6,0.4;...
                        0.6,0.6,0.9;...
                        0.9,0.6,0.6;...
                        0.9,0.9,0.6;...
                        0.9,0.9,0.9];

                    for iPDO = [3,5,9]
                        %
                        figure(iPDO)
                        switch iTime
                            case 1
                                Z = max(-15,log10(cSet{iPDO}));
                                contourf([0:size(cSet{iPDO},2)-1],[0:size(cSet{iPDO},1)-1],Z,30)
                                colorbar
                                hold on
                                scatter(True,DistSet{iPDO},100,'sk','filled')
                            case 2
                                scatter(True,DistSet{iPDO},100,'om','filled')
                            case 3
                                scatter(True,DistSet{iPDO},100,'^c','filled')
                        end
                        caxis([-15,0])
                        set(gca,'xlim',[0,500])
                        switch iPDO
                            case {1,2,3,10}
                                ylims = limsMcp;
                            case {4,5,6}
                                ylims = limsIfish;
                            case {7,8,9}
                                ylims = limsImcp;
                        end

                        xlabel('Number observed (FISH)')
                        ylabel(['Number observed (',labSet{iPDO},')'])
                        title(titles{iPDO})
                        
                        set(gca,'fontsize',16,'xlim',[0,500],'ylim',[0,ylims(2)])

                        %%
                        figure(1234);
                        Nmax = max(max(True),max(DistSet{iPDO}));
                        edges = -1:Nmax;
                        h1 = histogram(DistSet{iPDO},edges,'normalization','cdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(1,:)); hold on
                        cdfDistorted = h1.Values;
                        h3 = histogram(DistSet{iPDO},edges,'normalization','pdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(1,:)); hold on
                        pdfDistorted = h3.Values;

                        h2 = histogram(True,edges,'normalization','cdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(2,:));

                        figure(100*iTime+iPDO)
                        stairs([-1:Nmax-1],h2.Values,'LineWidth',4,'Color',col(2+3*(iTime-1),:)); hold on
                        stairs([-1:Nmax-1],h1.Values,'LineWidth',4,'Color',col(1+3*(iTime-1),:)); hold on
                        close(1234)

                        set(gca,'FontSize',16,'xlim',[-1,max(max(True),max(DistSet{iPDO}))],'ylim',[0,1.02])
                        xlabel('Number of mRNA')
                        ylabel('Cumulative Probability')
                        title(titles{iPDO})

                        figure(1000+100*iTime+iPDO)
                        edges = 0:binSize:Nmax;
                        histogram(True,edges,'normalization','pdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(2,:)); hold on
                        histogram(DistSet{iPDO},edges,'normalization','pdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(1,:)); hold on
                        set(gca,'FontSize',16,'xlim',[-1,max(max(True),max(DistSet{iPDO}))])
                        xlabel('Number of mRNA')
                        ylabel('Probability Mass')

                        figure(1234); clf
                        edges = 0:Nmax;
                        h1 = histogram(DistSet{iPDO},edges,'normalization','pdf');
                        h1 = h1.Values;
                        h2 = histogram(True,[0:size(cSet{iPDO},2)],'normalization','pdf');
                        h2 = h2.Values;
                        close(1234);

                        pDist = cSet{iPDO};
                        distPred = pDist*h2';

                        pDistBin=[];
                        for i = 1:floor(length(distPred)/binSize)
                            pDistBin(i) = sum(distPred((i-1)*binSize+1:i*binSize))/binSize;
                        end
                        figure(100*iTime+iPDO)
                        stairs([-1:length(distPred)],[0;cumsum(distPred);1],'--','LineWidth',4,'Color',col(3,:));

                        figure(1000+100*iTime+iPDO)
                        stairs([0:binSize:length(distPred)-binSize],(pDistBin),'LineWidth',4,'Color',col(3,:));
                        legend('FISH' ,labSet{iPDO},'FISH + PDO')
                        title(titles{iPDO})

                        try
                            results.KS(iPDO,iTime) = max(abs(cdfDistorted'-cumsum(distPred)));
                            results.EMD(iPDO,iTime) = sum(abs(cdfDistorted'-cumsum(distPred)));
                            tmp = pdfDistorted.*log(pdfDistorted./distPred);
                            results.KLD(iPDO,iTime) =  sum(tmp(isfinite(tmp)));
                            %                         results.KLD(iPDO,iTime) =
                        catch
                            nn=max(length(cdfDistorted),length(distPred));
                            cdfDistorted(end:nn+1) = cdfDistorted(end);
                            pdfDistorted(end+1:nn+1) = 0;
                            distPred(end+1:nn+1) = 0;
                            results.KS(iPDO,iTime) = max(abs(cdfDistorted'-cumsum(distPred)));
                            results.EMD(iPDO,iTime) = sum(abs(cdfDistorted'-cumsum(distPred)));
                            tmp = pdfDistorted.*log(pdfDistorted./distPred);
                            results.KLD(iPDO,iTime) =  sum(tmp(isfinite(tmp)));
                        end
                    end
                end
            end
        end

    case 32 % Fit model to data
        % For this fit, we will assume that the smFISH data is "true" and we will
        % fit once using the smFISH data and once with the MCP-GFP data.

        load(fileName,'ModelZero','pdoPars','chainResults')
        try
            if ~exist('ModelZero','var')
                error('missing ModelZero')
            end
        catch
            disp('Could not find previous model.  Starting from scratch.')
            [ModelZero,chainResults] = createDefaults(muLog10Prior,vars);
        end

        if isempty(vars)||vars.doFit==1
            parGuessesLocal =[];
            for iPDO = vars.fitCases
                parGuessesLocal = [parGuessesLocal;[ModelZero{iPDO}.parameters{vars.modelVarsToFit,2}]];
            end
            if ~isempty(vars.otherParGuesses)
                for iGuess = 1:length(vars.otherParGuesses)
                    try
                        modTmp = load(vars.otherParGuesses{iGuess},'ModelZero');
                        for iPDO = vars.fitCases
                            parGuessesLocal = [parGuessesLocal;[modTmp.ModelZero{iPDO}.parameters{vars.modelVarsToFit,2}]];
                        end
                    catch
                    end
                end
            end

            ModelZeroL = ModelZero(vars.fitCases);
            try
                chainResultsL = chainResults(vars.fitCases);
            catch
                for i = length(chainResults)+1:max(vars.fitCases)
                    chainResults{i} = [];
                end
                chainResultsL = chainResults(vars.fitCases);
            end
            errorDetected = zeros(1,length(ModelZeroL));
            parfor jPDO = 1:length(ModelZeroL)
                iPDO = vars.fitCases(jPDO);
                ModelZeroL{jPDO}.fspOptions.usePiecewiseFSP=false;
                ModelZeroL{jPDO}.fspOptions.initApproxSS=true;
                ModelZeroL{jPDO}.inputExpressions = {'ITrypt','(t<5)'};

                [ModelZeroL{jPDO},fitError{jPDO}] = runFittingFun(ModelZeroL{jPDO},vars,iPDO,...
                        pdoPars,muLog10Prior,sigLog10Prior,chainResultsL{jPDO},...
                        fitOptions,parGuessesLocal);
            end
            fitError

            if sum(errorDetected)~=0
                save([filename,'_ERROR'])
                error(['Error detected, data saved in ',filename,'_ERROR'])
            end

            ModelZero(vars.fitCases)=ModelZeroL;
            save(fileName,'ModelZero','-append')

        else

            parsTrue = ModelZero{1}.parameters;
            parsDistortedMCP = ModelZero{2}.parameters;
            parsCorrectedMCP = ModelZero{3}.parameters;
            parsCorrectedFISHIntens = ModelZero{4}.parameters;
            parsCorrectedMCPIntens = ModelZero{5}.parameters;
            T = table([parsTrue{:,2}]',[parsDistortedMCP{:,2}]',...
                [parsCorrectedMCP{:,2}]',[parsCorrectedFISHIntens{:,2}]',...
                [parsCorrectedMCPIntens{:,2}]');

            diffTrueDistortedMCP = exp(abs(log([parsDistortedMCP{:,2}]./[parsTrue{:,2}])))';
            diffTrueCorrectedMCP = exp(abs(log([parsCorrectedMCP{:,2}]./[parsTrue{:,2}])))';
            diffTrueCorrectedFISHIntens = exp(abs(log([parsCorrectedFISHIntens{:,2}]./[parsTrue{:,2}])))';
            diffTrueCorrectedMCPIntens = exp(abs(log([parsCorrectedMCPIntens{:,2}]./[parsTrue{:,2}])))';

            meanFoldDiffDistMCP = mean(diffTrueDistortedMCP(ModelZero{1}.fittingOptions.modelVarsToFit));
            meanFoldDiffCorectedMCP = mean(diffTrueCorrectedMCP(ModelZero{1}.fittingOptions.modelVarsToFit));
            meanFoldDiffCorectedFISHIntens = mean(diffTrueCorrectedFISHIntens(ModelZero{1}.fittingOptions.modelVarsToFit));
            meanFoldDiffCorectedMCPIntens = mean(diffTrueCorrectedMCPIntens(ModelZero{1}.fittingOptions.modelVarsToFit));

            for iPDO = vars.fitCases
                ModelZero{iPDO}.inputExpressions = {'ITrypt','(t<5)'};
                ModelZero{iPDO}.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
                ModelZero{iPDO}.fspOptions.fspTol = 0.0001;  % Set FSP error tolerance.
                ModelZero{iPDO}.tSpan = [-1e-6,0,5,18,300];

                ModelZero{iPDO}.fspOptions.bounds = vars.initialFspBounds;
                [FSPsoln,ModelZero{iPDO}.fspOptions.bounds] = ModelZero{iPDO}.solve;  % Solve the FSP analysis
                ModelZero{iPDO} = associateData(ModelZero{iPDO},vars,iPDO,pdoPars,FSPsoln,6000);
                ModelZero{iPDO}.tSpan = sort(unique([ModelZero{iPDO}.initialTime,5,ModelZero{iPDO}.dataSet.times]));

                [~,~,fitSolution{iPDO}] = ModelZero{iPDO}.computeLikelihood;

                results.fitLikelihoods(iPDO,:) = (fitSolution{iPDO}.DataLoadingAndFittingTabOutputs.V_LogLk-...
                    fitSolution{iPDO}.DataLoadingAndFittingTabOutputs.perfectMod);
                results.fitKS(iPDO,:) = fitSolution{iPDO}.DataLoadingAndFittingTabOutputs.V_KS;
            end

            for iPDO = 1:5
                ModelZero{iPDO}.makeFitPlot(fitSolution{iPDO},binSize,[iPDO,10+iPDO,103,104])  % Fit total
            end

            ModelTMP1 = ModelZero{1};
            for iPDO = 2:5
                ModelTMP1.parameters = ModelZero{iPDO}.parameters;
                ModelTMP1.fspOptions.fspTol = 0.0001;  % Set FSP error tolerance.
                [~,ModelTMP1.fspOptions.bounds] = ModelTMP1.solve;  % Solve the FSP analysis
                [~,~,fitSolution{4+iPDO}] = ModelTMP1.computeLikelihood;
                results.fitLikelihoods(4+iPDO,:) = (fitSolution{4+iPDO}.DataLoadingAndFittingTabOutputs.V_LogLk-...
                    fitSolution{4+iPDO}.DataLoadingAndFittingTabOutputs.perfectMod);
                results.fitKS(4+iPDO,:) = fitSolution{4+iPDO}.DataLoadingAndFittingTabOutputs.V_KS;

                ModelTMP1.makeFitPlot(fitSolution{4+iPDO},binSize,[iPDO,10+iPDO,103,104]);
            end
            results.fitSolution=fitSolution;
        end

    case 3 % Compute FIM for best parameter fit
        load(fileName,'ModelZero','pdoPars')
        ModelZeroL = ModelZero(vars.fitCases);
        parfor jPDO = 1:length(ModelZeroL)
            iPDO = vars.fitCases(jPDO);
            ModelZeroL{jPDO}.fspOptions.fspTol = 1e-6;  % Set FSP error tolerance.
            ModelZeroL{jPDO}.solutionScheme = 'FSP'; % Set solutions scheme to FSP
            ModelZeroL{jPDO}.fspOptions.bounds = vars.initialFspBounds;
            [FSPsoln,ModelZeroL{jPDO}.fspOptions.bounds] = ModelZeroL{jPDO}.solve;  % Solve the problem
            ModelZeroL{jPDO} = associateData(ModelZeroL{jPDO},vars,iPDO,pdoPars,FSPsoln);

            ModelZeroL{jPDO}.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
            ModelZeroL{jPDO}.sensOptions.solutionMethod = 'finiteDifference';

            ModelZeroL{jPDO}.tSpan = unique([[-1e-6,5],...
                ModelZeroL{jPDO}.dataSet.times(ModelZeroL{jPDO}.fittingOptions.timesToFit)]);

            [sensSoln,ModelZeroL{jPDO}.fspOptions.bounds] = ModelZeroL{jPDO}.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

            try
                fimResultsL{jPDO} = ModelZeroL{jPDO}.computeFIM(sensSoln.sens);
            catch
                if iPDO>=3
                    ModelZeroL{jPDO}.fspOptions.fspTol = 1e-7;  % Set FSP error tolerance.
                    ModelZeroL{jPDO}.solutionScheme = 'FSP'; % Set solutions scheme to FSP
                    [FSPsoln,ModelZeroL{jPDO}.fspOptions.bounds] = ModelZeroL{jPDO}.solve;  % Solve the problem
                    ModelZeroL{jPDO}.pdoOptions.PDO = ModelZeroL{jPDO}.generatePDO(ModelZeroL{jPDO}.pdoOptions,...
                        ModelZeroL{jPDO}.pdoOptions.props.PDOpars,FSPsoln.fsp,false);
                end
                fimResultsL{jPDO} = ModelZeroL{jPDO}.computeFIM(sensSoln.sens);
            end

            cellCounts = zeros(1,length(ModelZeroL{jPDO}.tSpan));
            for it = 1:length(ModelZeroL{jPDO}.dataSet.times)
                if ModelZeroL{jPDO}.fittingOptions.timesToFit(it)
                    [~,jm] = min(abs(ModelZeroL{jPDO}.tSpan-ModelZeroL{jPDO}.dataSet.times(it)));
                    cellCounts(jm) = ModelZeroL{jPDO}.dataSet.nCells(it);
                end
            end

            FIMZeroL{jPDO} = ModelZeroL{jPDO}.evaluateExperiment(fimResultsL{jPDO},cellCounts);

            FIMZeroLogL{jPDO} = diag([ModelZeroL{jPDO}.parameters{ModelZeroL{jPDO}.fittingOptions.modelVarsToFit,2}])*...
                FIMZeroL{jPDO}(ModelZeroL{jPDO}.fittingOptions.modelVarsToFit,ModelZeroL{jPDO}.fittingOptions.modelVarsToFit)*...
                diag([ModelZeroL{jPDO}.parameters{ModelZeroL{jPDO}.fittingOptions.modelVarsToFit,2}]);
            covLogL{jPDO} = FIMZeroLogL{jPDO}^-1;
        end
        FIMZero(vars.fitCases)=FIMZeroL;
        FIMZeroLog(vars.fitCases)=FIMZeroLogL;
        fimResults(vars.fitCases)=fimResultsL;
        covLog(vars.fitCases)=covLogL;
        save(fileName,'FIMZero','FIMZeroLog','fimResults','covLog','-append')

    case 4 % Run Met Hast to quantify uncertainty on initial parameters
        if isempty(vars)||vars.doFit==1
            load(fileName,'ModelZero','FIMZeroLog','pdoPars')
            for iPDO = vars.fitCases
                indsPars = ModelZero{iPDO}.fittingOptions.modelVarsToFit;
                ModelZero{iPDO}.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior(indsPars)).^2./(2*sigLog10Prior(indsPars).^2);
            end
            ModelZeroL = ModelZero(vars.fitCases);
            FIMZeroLogL = FIMZeroLog(vars.fitCases);
            parfor jPDO = 1:length(ModelZeroL)
                iPDO = vars.fitCases(jPDO);

                fimLog = FIMZeroLogL{jPDO} + diag((1./vars.sigLog10Prior(indsPars)).^2/log(10)^2);
                covLog = (fimLog)^-1;

                MHOptions = struct('numberOfSamples',vars.nMH,'burnIn',0,'thin',1,'saveFile',['TMPmh_',num2str(iPDO),'_',fileName],...
                    'useFIMforMetHast',false);
                if length(vars.mhScaling)==1
                    proposalWidthScale = vars.mhScaling
                else
                    proposalWidthScale = vars.mhScaling(iPDO)
                end
                MHOptions.proposalDistribution  = @(x)mvnrnd(x,proposalWidthScale*(covLog+covLog')/2);

                % Solve initial model and compute likelihood
                ModelZeroL{jPDO}.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
                ModelZeroL{jPDO}.fspOptions.fspTol = 0.00001;  % Set FSP error tolerance.
                ModelZeroL{jPDO}.fspOptions.bounds = vars.initialFspBounds;
                [FSPsoln,ModelZeroL{jPDO}.fspOptions.bounds] = ModelZeroL{jPDO}.solve;  % Solve the problem
                ModelZeroL{jPDO} = associateData(ModelZeroL{jPDO},vars,iPDO,pdoPars,FSPsoln);
                ModelZeroL{jPDO}.tSpan = sort(unique([ModelZeroL{jPDO}.initialTime,5,ModelZeroL{jPDO}.dataSet.times(ModelZeroL{jPDO}.fittingOptions.timesToFit)]));
                
                % set weak tolerance for search.
                ModelZeroL{jPDO}.fspOptions.fspTol = inf;  % Set FSP error tolerance.
                mlikelihood = -ModelZeroL{jPDO}.computeLikelihood([],FSPsoln.stateSpace);

                % run MH starting at current parameter set
                parGuess = [ModelZeroL{jPDO}.parameters{ModelZeroL{jPDO}.fittingOptions.modelVarsToFit,2}];
                [~,~,chainResultsL{jPDO}] = ModelZeroL{jPDO}.maximizeLikelihood(parGuess',MHOptions,'MetropolisHastings');

            end
            
            chainResults(vars.fitCases)=chainResultsL;
            save(fileName,'chainResults','-append')
            for i=vars.fitCases
                disp('MH Completed successfully. Deleting temporary files.')
                delete(['TMPmh_',num2str(i),'_',fileName,'.mat'])
            end

        elseif vars.doFit==0
            load(fileName,'chainResults','covLog','ModelZero')
            for iPDO = vars.fitCases
                log10MLEPars = log10([ModelZero{iPDO}.parameters{ModelZero{iPDO}.fittingOptions.modelVarsToFit,2}]);
                smplDone = chainResults{iPDO}.mhSamples(:,:);
                valDone = chainResults{iPDO}.mhValue;
                smplDone = smplDone(valDone~=0,:);
                valDone = valDone(valDone~=0);

                [valDoneSorted,J] = sort(valDone);
                smplDone = smplDone(J,:);

                covPRIOR = diag(sigLog10Prior(ModelZero{iPDO}.fittingOptions.modelVarsToFit));
                muPRIOR = muLog10Prior(ModelZero{iPDO}.fittingOptions.modelVarsToFit);
                sigPRIOR = sigLog10Prior(ModelZero{iPDO}.fittingOptions.modelVarsToFit);
                covFIM = covLog{iPDO}/(log(10)^2);

                if vars.covWithPrior
                    % Add effect of prior to FIM.
                    FIMPrior = diag(1./sigPRIOR.^2);% + logDeltPars*logDeltPars';
                    covFIM = (covFIM^-1+FIMPrior)^-1;
                end

                covMLE = cov(smplDone/log(10));
                stdMLE = sqrt(diag(covMLE))';
                stdFIM = sqrt(diag(covFIM))';
                stdsFIM(iPDO,:) = stdFIM;
                stdsMLE(iPDO,:) = stdMLE;
                detsCovFIM(iPDO) = det(covFIM);
                detsCovMLE(iPDO) = det(covMLE);

                [evecsFIMinv,evalsFIMinv] = eig(covFIM);
                [evecsMH,evalsMH] = eig(covMLE);
                

                figure(123);
                subplot(1,5,iPDO);

                if nargout>=1
                    results.covMLE{iPDO} = covMLE;
                    results.detCovMLE = detsCovMLE;
                    results.muMLE(iPDO,:) = mean(exp(smplDone));
                    results.MLE(iPDO,:) = log10MLEPars;
                    results.sigMLE(iPDO,:) = std(exp(smplDone));
                    results.evecsFIMinv{iPDO} = evecsFIMinv;
                    results.evalsFIMinv{iPDO} = evalsFIMinv;
                    results.evecsMH{iPDO} = evecsMH;
                    results.evalsMH{iPDO} = evalsMH;
                else

                    for ipar = 1:size(covFIM,1)
                        errorbar(muPRIOR(ipar),ipar+0.15,covPRIOR(ipar,ipar),'cs','horizontal','LineWidth',3)
                        hold on
                        errorbar(log10MLEPars(ipar),ipar-0.15,stdMLE(ipar),'ms','horizontal','LineWidth',3)
                        errorbar(log10MLEPars(ipar),ipar,sqrt(covFIM(ipar,ipar)),'ks','horizontal','LineWidth',3)
                        if iPDO==1
                            parsTotal = log10MLEPars;
                        end
                        plot(parsTotal(ipar)*[1,1],ipar+[-0.3,0.3],'k--')
                    end
                    set(gca,'ylim',[0,6])

                    figure(1+(iPDO-1)*3)
                    plot(valDone)
                    figure(2+(iPDO-1)*3)

                    for i =1:size(smplDone,2)-1
                        for j = i+1:size(smplDone,2)
                            %% Plot of Log10 MH results
                            figure(2+(iPDO-1)*3)
                            subplot(size(smplDone,2)-1,size(smplDone,2)-1,(i-1)*(size(smplDone,2)-1)+j-1)
                            if vars.mleScatterPlots
                                scatter(smplDone(:,j)/log(10),smplDone(:,i)/log(10),6,valDoneSorted,'filled'); hold on;
                                plot(smplDone(end,j)/log(10),smplDone(end,i)/log(10),'bs','markerfacecolor','b');
                            end
                            par0 = mean(smplDone(:,[j,i])/log(10));

                            covLogI = covFIM([j,i],[j,i]);

                            ssit.parest.ellipse(log10MLEPars([j,i]),icdf('chi2',0.9,2)*covLogI,'k-','linewidth',3)
                            hold on
                            clear cov
                            cov12 = cov(smplDone(:,j)/log(10),smplDone(:,i)/log(10));
                            ssit.parest.ellipse(par0,icdf('chi2',0.9,2)*cov12,'m--','linewidth',2)

                            ssit.parest.ellipse(muPRIOR([j,i]),icdf('chi2',0.9,2)*covPRIOR([j,i],[j,i]),'c-','linewidth',3)
                            plot(log10MLEPars(j),log10MLEPars(i),'kx','MarkerSize',8)

                            if vars.mleScatterPlots
                                figure(3+(iPDO-1)*3)
                                subplot(size(smplDone,2)-1,size(smplDone,2)-1,(i-1)*(size(smplDone,2)-1)+j-1)
                                plot(exp(smplDone(:,j)),exp(smplDone(:,i)),'ro')
                            end

                        end
                    end
                end
            end
        end
    case 5 % Sample from MH Results and Compute FIM at all time points
        clear ModelOne sensSolnOne fimResultsOne
        load(fileName,'ModelZero','chainResults','pdoPars')
        nSamples = vars.nFIMsamples;

        chainResultsL = chainResults(vars.fitCases);
        ModelZeroL = ModelZero(vars.fitCases);
        parfor jPDO = 1:length(vars.fitCases)
            iPDO = vars.fitCases(jPDO);
            valDone = chainResults{jPDO}.mhValue;
            smplDone = chainResultsL{jPDO}.mhSamples(:,:);
            smplDone = smplDone(valDone~=0,:);
            ismpls = floor(linspace(1,size(smplDone,1),nSamples));
            parsetsL{jPDO} = repmat([ModelZeroL{jPDO}.parameters{:,2}],nSamples,1);
            parsetsL{jPDO}(:,vars.modelVarsToFit) = exp(smplDone(ismpls,:));
            parsetsL{jPDO}(1,:) = [ModelZeroL{jPDO}.parameters{:,2}];

            for i=1:nSamples
                [i,jPDO]

                ModelOne = ModelZeroL{jPDO};
                ModelOne.tSpan = sort(unique([-1e-6,0:6:1200,5]));
                ModelOne.parameters(:,2) = num2cell(parsetsL{jPDO}(i,:));

                ModelOne.solutionScheme = 'FSP'; % Set solutions scheme to FSP Sensitivity
                ModelOne.fspOptions.fspTol = 1e-7;
                [FSPsoln,ModelOne.fspOptions.bounds] = ModelOne.solve;  % Solve the FSP analysis
                ModelOne = associateData(ModelOne,vars,iPDO,pdoPars,FSPsoln);

                ModelOne.fspOptions.fspTol = inf;
                ModelOne.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
                ModelOne.sensOptions.solutionMethod = 'finiteDifference';
                
                ModelOne.tSpan = sort(unique([-1e-6,0:6:1200,5]));
                sensSolnOne = ModelOne.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

                try
                    fimResultsOneL{i,jPDO} = ModelOne.computeFIM(sensSolnOne.sens);
                catch
                    if iPDO>=3
                        ModelOne.fspOptions.fspTol = 1e-7;  % Set FSP error tolerance.
                        ModelOne.solutionScheme = 'FSP'; % Set solutions scheme to FSP
                        [FSPsoln,ModelOne.fspOptions.bounds] = ModelOne.solve;  % Solve the problem
                        ModelOne.pdoOptions.PDO = ModelOne.generatePDO(ModelOne.pdoOptions,...
                            ModelOne.pdoOptions.props.PDOpars,FSPsoln.fsp,false);
                    end
                    fimResultsOneL{i,jPDO} = ModelOne.computeFIM(sensSolnOne.sens);
                end
            end
        end
        fimResultsOne(:,vars.fitCases) = fimResultsOneL;
        parsets(vars.fitCases) = parsetsL;

        disp('Saving FIM v time results');
        save(fileName,'fimResultsOne','parsets','-append')

    case 6 % Use FIM to Estimate Experiment Designs.
        load(fileName,'ModelZero','fimResultsOne','parsets')

        designChoice = 'DetCovDesign';
        minormax = 1;
        titltxt = 'det(FIM^{-1})';

        DetCovDesign=[];
        DetDesign=[];
        EDesign=[];
        TDesign=[];
        gDesign=[];

        for jpar = size(fimResultsOne,1):-1:1
            for iPDO = vars.fitCases
                for it=size(fimResultsOne{1,1},1):-1:1
                    fimResultsRed{jpar,iPDO}{it,1} = ...
                        diag(parsets{iPDO}(jpar,vars.modelVarsToFit))*...
                        fimResultsOne{jpar,iPDO}{it}(vars.modelVarsToFit,vars.modelVarsToFit)*...
                        diag(parsets{iPDO}(jpar,vars.modelVarsToFit));
                end
            end
        end

        for jpar = 1:size(fimResultsRed,1)
            for iPDO = vars.fitCases
                for it = 2:size(fimResultsRed{1,1},1)
                    nCellsDesign = zeros(1,size(fimResultsRed{1,1},1));
                    nCellsDesign(2) = vars.FIMnCells(1);
                    nCellsDesign(6) = vars.FIMnCells(2);
                    nCellsDesign(53) = vars.FIMnCells(3);
                    nCellsDesign(it) = nCellsDesign(it)+100;
                    [fimAll] = ModelZero{1}.evaluateExperiment(fimResultsRed{jpar,iPDO},nCellsDesign);
                    
                    if vars.covWithPrior
                        fimAll = fimAll + diag((1./vars.sigLog10Prior(vars.modelVarsToFit)).^2/log(10)^2);
                    end

                    % Change base of FIM to log10.
                    fimAll = fimAll*log(10)^2;

                    covAll = inv(fimAll);

                    DetDesign(jpar,iPDO,it) = det(fimAll);
                    DetCovDesign(jpar,iPDO,it) = det(covAll);
                    EDesign(jpar,iPDO,it) = min(eig(fimAll));
                    TDesign(jpar,iPDO,it) = trace(fimAll);
                    if ~isnan(covAll)
                        gDesign(jpar,iPDO,it) = covAll(end,end);
                    else
                        gDesign(jpar,iPDO,it) = NaN;
                    end
                end
            end
        end
        %%

        eval(['design = ',designChoice,';']);

        logFIMCritIDeal = log10(squeeze(design(:,1,2:end)));
        logmnIdeal = mean(logFIMCritIDeal);
        logstdIdeal = std(logFIMCritIDeal);
 
        logFIMCritPdoFishInten = log10(squeeze(design(:,4,2:end)));
        logmnPdoFishInten = mean(logFIMCritPdoFishInten);
        logstdPdoFishInten = std(logFIMCritPdoFishInten);

        logFIMCritPdoMcpSpots = log10(squeeze(design(:,3,2:end)));
        logmnPdoMcpSpots= mean(logFIMCritPdoMcpSpots);
        logstdPdoMcpSpots = std(logFIMCritPdoMcpSpots);

        logFIMCritPdoMcpInten = log10(squeeze(design(:,5,2:end)));
        logmnPdoMcpInten= mean(logFIMCritPdoMcpInten);
        logstdPdoMcpInten = std(logFIMCritPdoMcpInten);

        if length(logFIMCritIDeal)==202
            xvals = sort(unique([1,6:6:1200,5]));
        else
            xvals = sort(unique([-1e-6,0:6:1200,5]));
        end

        figure

        if ~isfield(vars,'FIMShowErrorBars')||vars.FIMShowErrorBars
            makePatch(log10(xvals),logmnIdeal,logstdIdeal,[0.8,0.4,0.8]); hold on
            makePatch(log10(xvals),logmnPdoMcpSpots,logstdPdoMcpSpots,[0.4,0.4,0.4])
            makePatch(log10(xvals),logmnPdoFishInten,logstdPdoFishInten,[0.8,0.8,0.4])
            makePatch(log10(xvals),logmnPdoMcpInten,logstdPdoMcpInten,[0.8,0.4,0.4])
            plot(log10(xvals),logmnIdeal,...
                'color','k','LineWidth',3); hold on
            plot(log10(xvals),logmnPdoMcpSpots,...
                'color','k','LineWidth',3); hold on
            plot(log10(xvals),logmnPdoFishInten,...
                'color','k','LineWidth',3); hold on
            plot(log10(xvals),logmnPdoMcpInten,...
                'color','k','LineWidth',3); hold on
        else
            plot(xvals,logmnIdeal,...
                'color',[0.8,0.4,0.8],'LineWidth',3); hold on
            plot(xvals,logmnPdoMcpSpots,...
                'color',[0.4,0.8,0.4],'LineWidth',3);
            plot(xvals,logmnPdoFishInten,...
                'color',[0.8,0.8,0.4],'LineWidth',3);
            plot(xvals,logmnPdoMcpInten,...
                'color',[0.8,0.4,0.4],'LineWidth',3);
        end


        ylim = get(gca,'ylim');
        tt = xvals;
        if minormax==1
            [~,j1] = min(logmnIdeal);
            [~,j2] = min(logmnPdoFishInten);
            [~,j3] = min(logmnPdoMcpSpots);
            [~,j4] = min(logmnPdoMcpInten);
        else
            [~,j1] = max(logmnIdeal);
            [~,j2] = max(logmnPdoFishInten);
            [~,j3] = max(logmnPdoMcpSpots);
            [~,j4] = max(logmnPdoMcpInten);
        end
        tt([j1,j3,j2,j4])
        plot(log10(tt(j1))*[1,1],ylim+[-10,10],'color',[0.8,0.4,0.8],'LineWidth',3)
        plot(log10(tt(j3))*[1,1],ylim+[-10,10],'color',[0.4,0.8,0.4],'LineWidth',3)
        plot(log10(tt(j2))*[1,1],ylim+[-10,10],'color',[0.8,0.8,0.4],'LineWidth',3)
        plot(log10(tt(j4))*[1,1],ylim+[-10,10],'color',[0.8,0.4,0.4],'LineWidth',3)

        plot(log10(tt(j1))*[1,1],ylim+[-10,10],'k--','LineWidth',3)
        plot(log10(tt(j3))*[1,1],ylim+[-10,10],'k--','LineWidth',3)
        plot(log10(tt(j2))*[1,1],ylim+[-10,10],'k--','LineWidth',3)
        plot(log10(tt(j4))*[1,1],ylim+[-10,10],'k--','LineWidth',3)

        set(gca,'FontSize',16); xlabel('Time (min)'); ylabel(titltxt)
        legend('Ideal','MCP-Spots+PDO','FISH-Intens+PDO','MCP-Intens+PDO')

end
cleanUpFile(fileName)
end

%%  Functions
function [logL,P] = findError(lambda,True,Distorted,lims)
arguments
    lambda
    True
    Distorted
    lims
end
% Computes likelihood of observed data given the model of affine poisson
% extra spot counting and probability of measurmeent failure.
logL = sum(max(-400,log(pdf('poiss',Distorted,max(lambda(1),lambda(2)+lambda(3)*True)))));
logL = logL-1e4*(lambda(1)<0);

if nargout>=2
%     Nmax = max(max(max(True,Distorted)));
%     Np = max(1,lambda(2)+lambda(3)*Nmax);
%     Np = max(Nmax,ceil(Np+5*sqrt(Np)));
    Np = lims(2);
    maxx = 100*ceil(max(True)/100);
    P = zeros(Np+1,maxx+1);
    for xi = 0:maxx
        % Affine Poisson Gain followed by binomial loss
        P(1:Np+1,xi+1) = pdf('poiss',[0:Np]',max(lambda(1),lambda(2)+lambda(3)*xi));
    end
end
end

function [logL,P] = findErrorFail(lambda,True,Distorted,lims)
arguments
    lambda
    True
    Distorted
    lims
end
% Computes likelihood of observed data given the model of affine poisson
% extra spot counting and probability of measurmeent failure.
% This PDO also allows for the measurement to fail (i.e., return zero spots)
% with probability lambda(4).
likelikood = lambda(4)*(Distorted==0)+...
    (1-lambda(4))*pdf('poiss',Distorted,max(lambda(1),lambda(2)+lambda(3)*True));

logL = sum(max(-400,log(likelikood)));
logL = logL-1e4*(lambda(1)<0)-1e4*(lambda(4)<0)-1e4*(lambda(4)>1);

if nargout>=2
%     Nmax = max(max(max(True,Distorted)));
%     Np = max(1,lambda(2)+lambda(3)*Nmax);
%     Np = max(Nmax,ceil(Np+5*sqrt(Np)));
    Np = lims(2);
    maxx = 100*ceil(max(True)/100);
    P = zeros(Np+1,maxx+1);
    for xi = 0:maxx
        % Affine Poisson Gain followed by binomial loss
        P(1:Np+1,xi+1) = pdf('poiss',[0:Np]',max(lambda(1),lambda(2)+lambda(3)*xi));
    end
    P = (1-lambda(4))*P;
    P(1,:) = P(1,:)+lambda(4);
end
end

function [logL,P,edges] = findErrorNorm(lambda,True,Distorted,lims)
arguments
    lambda
    True
    Distorted
    lims
end
% Computes likelihood of observed data given the PDO model of intensity,
% where there is a gaussian background and a gausian intensity per spot.
dif = Distorted - (lambda(1)+lambda(3)*True);
sig2 = lambda(2)+lambda(4)*True;
% logpdf = -1/2*log(2*pi*sig2)-dif.^2./(2*sig2);
% I = llogpdf>-4;
logcdf = log(cdf("Normal",dif+.5,0,sqrt(sig2))-cdf("Normal",dif-.5,0,sqrt(sig2)));
logL = sum(max(-400,min(0,logcdf)));
logL = logL-1e4*((lambda(1)<0)+(lambda(2)<0)+(lambda(3)<0)+(lambda(4)<0));

if nargout>=2
    nIntensBins = lims(2)+1;
    edges = [0:nIntensBins];
    maxx = 100*ceil(max(True)/100);
    P = zeros(nIntensBins+1,maxx);
    for i=0:maxx
        mn = (lambda(1)+lambda(3)*i);
        sig = sqrt(lambda(2)+lambda(4)*i);
%         P(:,i+1) = pdf('norm',edges,mn,sig);
        P(:,i+1) = cdf('norm',edges+.5,mn,sig)-cdf('norm',edges-.5,mn,sig);
        if sum(P(:,i+1)) > 1
           P(:,i+1)=P(:,i+1)/sum(P(:,i+1));
        end
    end
end
end

function [logL,P,edges] = findErrorNormFail(lambda,True,Distorted,lims)
arguments
    lambda
    True
    Distorted
    lims
end
% Computes likelihood of observed data given the PDO model of intensity,
% where there is a gaussian background and a gausian intensity per spot.
dif = Distorted - (lambda(1)+lambda(3)*True);
sig2 = lambda(2)+lambda(4)*True;
difFail = Distorted - (lambda(1));
sig2Fail = lambda(2);

% logpdf = -1/2*log(2*pi*sig2)-dif.^2./(2*sig2);
% I = llogpdf>-4;
probs = cdf("Normal",dif+.5,0,sqrt(sig2))-cdf("Normal",dif-.5,0,sqrt(sig2));
probsFail = cdf("Normal",difFail+.5,0,sqrt(sig2Fail))-cdf("Normal",difFail-.5,0,sqrt(sig2Fail));

probs = lambda(5)*probsFail+(1-lambda(5))*probs;
logcdf = log(probs);
logL = sum(max(-400,min(0,logcdf)));
logL = logL-1e4*((lambda(1)<0)+(lambda(2)<0)+(lambda(3)<0)+(lambda(4)<0)...
    +(lambda(5)<0)+(lambda(5)>1));

if nargout>=2
    nIntensBins = lims(2)+1;
    edges = [0:nIntensBins];
    maxx = 100*ceil(max(True)/100);
    P = zeros(nIntensBins+1,maxx);
    for i=0:maxx
        mn = (lambda(1)+lambda(3)*i);
        sig = sqrt(lambda(2)+lambda(4)*i);
        mnFail = lambda(1);
        sigFail = sqrt(lambda(2));
        %         P(:,i+1) = pdf('norm',edges,mn,sig);
        P(:,i+1) = (1-lambda(5))*(cdf('norm',edges+.5,mn,sig)-cdf('norm',edges-.5,mn,sig))+...
            lambda(5)*(cdf('norm',edges+.5,mnFail,sigFail)-cdf('norm',edges-.5,mnFail,sigFail));
        if sum(P(:,i+1)) > 1
            P(:,i+1)=P(:,i+1)/sum(P(:,i+1));
        end
    end
end
end

function [logL,P,edges] = findErrorNormSpurrious(lambda,True,Distorted,lims)
arguments
    lambda
    True
    Distorted
    lims
end
% Computes likelihood of observed data given the PDO model of intensity,
% where there is a gaussian background and a gausian intensity per spot.
dif = Distorted - (lambda(1)+lambda(3)*True);
sig2 = lambda(2)+lambda(4)*True;
% difFail = Distorted - lambda(6) + lambda(3)*True;
% sig2Fail = lambda(7)+lambda(2)+lambda(4)*True;
difFail = Distorted - lambda(6) + 0*True;
sig2Fail = lambda(7) + 0*True;

% logpdf = -1/2*log(2*pi*sig2)-dif.^2./(2*sig2);
% I = llogpdf>-4;
probs = cdf("Normal",dif+.5,0,sqrt(sig2))-cdf("Normal",dif-.5,0,sqrt(sig2));
probs(abs(dif)>3*sqrt(sig2)) = pdf("Normal",dif(abs(dif)>3*sqrt(sig2)),0,sqrt(sig2(abs(dif)>3*sqrt(sig2))));

probsFail = cdf("Normal",difFail+.5,0,sqrt(sig2Fail))-cdf("Normal",difFail-.5,0,sqrt(sig2Fail));
probsFail(abs(difFail)>3*sqrt(sig2Fail)) = ...
    pdf("Normal",difFail(abs(difFail)>3*sqrt(sig2Fail)),0,sqrt(sig2Fail(abs(difFail)>3*sqrt(sig2Fail))));

% Project negative values to zero bin
probs(Distorted==0) = cdf("Normal",dif(Distorted==0)+.5,0,sqrt(sig2(Distorted==0)));
probsFail(Distorted==0) = cdf("Normal",difFail(Distorted==0)+.5,0,sqrt(sig2Fail(Distorted==0)));

NFail = 10;
probs(True<=NFail) = lambda(5)*probsFail(True<=NFail)+(1-lambda(5))*probs(True<=NFail);

logcdf = log(probs);
logL = sum(max(-400,min(0,logcdf)));

if (lambda(1)<0)||(lambda(2)<0)||(lambda(3)<0)...
        ||(lambda(4)<0)||(lambda(5)<0)||(lambda(5)>1)||(lambda(7)<0)
    logL = -1e20;
end

if nargout>=2
    nIntensBins = lims(2)+1;
    edges = [0:nIntensBins];
    maxx = 100*ceil(max(True)/100);
    P = zeros(nIntensBins+1,maxx);
    for i=0:maxx
        mn = (lambda(1)+lambda(3)*i);
        sig = sqrt(lambda(2)+lambda(4)*i);
%         mnFail = lambda(6) + lambda(3)*i;
%         sigFail = sqrt(lambda(7)+lambda(2)+lambda(4)*i);
        mnFail = lambda(6);
        sigFail = sqrt(lambda(7));
        %         P(:,i+1) = pdf('norm',edges,mn,sig);
        if i<=NFail
            P(:,i+1) = (1-lambda(5))*(cdf('norm',edges+.5,mn,sig)-cdf('norm',edges-.5,mn,sig))+...
                lambda(5)*(cdf('norm',edges+.5,mnFail,sigFail)-cdf('norm',edges-.5,mnFail,sigFail));
            P(10:end,i+1) = (1-lambda(5))*pdf('norm',edges(10:end),mn,sig)+...
                lambda(5)*pdf('norm',edges(10:end),mnFail,sigFail);
            P(1,i+1) = (1-lambda(5))*(cdf('norm',.5,mn,sig))+...
                lambda(5)*(cdf('norm',.5,mnFail,sigFail));
        else
            P(1:9,i+1) = cdf('norm',edges(1:9)+.5,mn,sig)-cdf('norm',edges(1:9)-.5,mn,sig);
            P(10:end,i+1) = pdf('norm',edges(10:end),mn,sig);
            P(1,i+1) = cdf('norm',.5,mn,sig);
        end
%         if sum(P(:,i+1)) > 1
%             P(:,i+1)=P(:,i+1)/sum(P(:,i+1));
%         end
    end
end
end

function [logL,P,edges] = findErrorNormAdd(lambda,True,Distorted,lims)
arguments
    lambda
    True
    Distorted
    lims
end
% Computes likelihood of observed data given the PDO model of intensity,
% where there is a gaussian background and a gausian intensity per spot.
dif = Distorted - True;
sig2 = lambda(2)+0*dif;
difFail = Distorted - True;
sig2Fail = lambda(3)+0*difFail;

probs = cdf("Normal",dif+.5,0,sqrt(sig2))-cdf("Normal",dif-.5,0,sqrt(sig2));
probs(abs(dif)>3*sqrt(sig2)) = pdf("Normal",dif(abs(dif)>3*sqrt(sig2)),0,sqrt(sig2(abs(dif)>3*sqrt(sig2))));

probsFail = cdf("Normal",difFail+.5,0,sqrt(sig2Fail))-cdf("Normal",difFail-.5,0,sqrt(sig2Fail));
probsFail(abs(difFail)>3*sqrt(sig2Fail)) = ...
    pdf("Normal",difFail(abs(difFail)>3*sqrt(sig2Fail)),0,sqrt(sig2Fail(abs(difFail)>3*sqrt(sig2Fail))));

% Project negative values to zero bin
probs(Distorted==0) = cdf("Normal",dif(Distorted==0)+.5,0,sqrt(sig2(Distorted==0)));
probsFail(Distorted==0) = cdf("Normal",difFail(Distorted==0)+.5,0,sqrt(sig2Fail(Distorted==0)));

probs = lambda(1)*probs+(1-lambda(1))*probsFail;

logcdf = log(probs);
logL = sum(max(-400,min(0,logcdf)));

if (lambda(1)<0)||(lambda(1)>1)||(lambda(2)<0)...
        ||(lambda(3)<0)
    logL = -1e20;
end

if nargout>=2
    nIntensBins = lims(2)+1;
    edges = [0:nIntensBins];
    maxx = 100*ceil(max(True)/100);
    P = zeros(nIntensBins+1,maxx);
    for i=0:maxx
        mn = i;
        sig = sqrt(lambda(2));
%         mnFail = lambda(6) + lambda(3)*i;
%         sigFail = sqrt(lambda(7)+lambda(2)+lambda(4)*i);
        mnFail = i;
        sigFail = sqrt(lambda(3));
        %         P(:,i+1) = pdf('norm',edges,mn,sig);
%         if i<=NFail
            P(:,i+1) = lambda(1)*(cdf('norm',edges+.5,mn,sig)-cdf('norm',edges-.5,mn,sig))+...
                (1-lambda(1))*(cdf('norm',edges+.5,mnFail,sigFail)-cdf('norm',edges-.5,mnFail,sigFail));
            P(10:end,i+1) = (lambda(1))*pdf('norm',edges(10:end),mn,sig)+...
                (1-lambda(1))*pdf('norm',edges(10:end),mnFail,sigFail);
            P(1,i+1) = lambda(1)*(cdf('norm',.5,mn,sig))+...
                (1-lambda(1))*(cdf('norm',.5,mnFail,sigFail));
%         else
%             P(1:9,i+1) = cdf('norm',edges(1:9)+.5,mn,sig)-cdf('norm',edges(1:9)-.5,mn,sig);
%             P(10:end,i+1) = pdf('norm',edges(10:end),mn,sig);
%             P(1,i+1) = cdf('norm',.5,mn,sig);
%         end
%         if sum(P(:,i+1)) > 1
%             P(:,i+1)=P(:,i+1)/sum(P(:,i+1));
%         end
    end
end
end

function ModelZero = associateData(ModelZero,vars,iPDO,pdoPars,FSPsoln,outputMax)
arguments
ModelZero
vars
iPDO
pdoPars
FSPsoln
outputMax =[];
end
ModelZero.pdoOptions= struct('unobservedSpecies',[],'PDO',[]); % Options for FIM analyses
obsvSpecies = find(strcmp(ModelZero.species,'x4'));
switch iPDO
    case 1 % FISH spots, no distortion model --- TRUE DATA
        ModelZero.pdoOptions.unobservedSpecies =  ModelZero.species(~contains(ModelZero.species,'x4'))';
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','number_spots_type_1'});
    case 2 % MCP-GFP Spot data, no distortion model
        ModelZero.pdoOptions.unobservedSpecies =  ModelZero.species(~contains(ModelZero.species,'x4'))';
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','number_spots_type_0'});
    case 3 % MCP-GFP Spot data, Gaussian + Spurrious PDO
        ModelZero.pdoOptions.unobservedSpecies =  [];
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','number_spots_type_0'});
        LL = pdoPars.lNormSpurriousMCPspots;
        ModelZero.pdoOptions.type = 'GaussSpurrious';
        ModelZero.pdoOptions.props.PDOpars = zeros(1,7*length(ModelZero.species));
        ModelZero.pdoOptions.props.PDOpars((obsvSpecies-1)*7+[1:7]) = LL;
        ModelZero.pdoOptions.props.pdoOutputRange = 3*ones(1,length(ModelZero.species));
        if isempty(outputMax)
            ModelZero.pdoOptions.props.pdoOutputRange(obsvSpecies) = ModelZero.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor.size(2);
        else
            ModelZero.pdoOptions.props.pdoOutputRange(obsvSpecies) = outputMax;
        end
        ModelZero.pdoOptions.PDO = ModelZero.generatePDO(ModelZero.pdoOptions,...
            ModelZero.pdoOptions.props.PDOpars,FSPsoln.fsp,false);
    case 4 % FISH Intensity Data data, Gaussian + Spurrious PDO
        ModelZero.pdoOptions.unobservedSpecies =  [];
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','nucIntens3'});
        LL = pdoPars.lNormSpurriousFISHintens;
        ModelZero.pdoOptions.type = 'GaussSpurrious';
        ModelZero.pdoOptions.props.PDOpars = zeros(1,7*length(ModelZero.species));
        ModelZero.pdoOptions.props.PDOpars((obsvSpecies-1)*7+[1:7]) = LL;
        ModelZero.pdoOptions.props.pdoOutputRange = 3*ones(1,length(ModelZero.species));
        if isempty(outputMax)
            ModelZero.pdoOptions.props.pdoOutputRange(obsvSpecies) = ModelZero.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor.size(2);
        else
            ModelZero.pdoOptions.props.pdoOutputRange(obsvSpecies) = outputMax;
        end
        ModelZero.pdoOptions.PDO = ModelZero.generatePDO(ModelZero.pdoOptions,...
            ModelZero.pdoOptions.props.PDOpars,FSPsoln.fsp,false);

    case 6 % MCP-GFP INTENSITY data, Gaussian
        ModelZero.pdoOptions.unobservedSpecies =  [];
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','nucIntens1'});
        LL = pdoPars.lNormTrue2MCPintens;
        ModelZero.pdoOptions.type = 'DiscretizedNormal';
        ModelZero.pdoOptions.props.PDOpars = zeros(1,4*length(ModelZero.species));
        ModelZero.pdoOptions.props.PDOpars((obsvSpecies-1)*4+[1:4]) = LL;
        ModelZero.pdoOptions.props.pdoOutputRange = 3*ones(1,length(ModelZero.species));
        if isempty(outputMax)
            ModelZero.pdoOptions.props.pdoOutputRange(obsvSpecies) = ModelZero.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor.size(2);
        else
            ModelZero.pdoOptions.props.pdoOutputRange(obsvSpecies) = outputMax;
        end
        ModelZero.pdoOptions.PDO = ModelZero.generatePDO(ModelZero.pdoOptions,...
            ModelZero.pdoOptions.props.PDOpars,FSPsoln.fsp,false);

    case 5 % MCP-GFP INTENSITY data, Gaussian + Spurrious PDO
        ModelZero.pdoOptions.unobservedSpecies =  [];
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','nucIntens1'});
        LL = pdoPars.lNormSpurriousMCPintens;
        ModelZero.pdoOptions.type = 'GaussSpurrious';
        ModelZero.pdoOptions.props.PDOpars = zeros(1,7*length(ModelZero.species));
        ModelZero.pdoOptions.props.PDOpars((obsvSpecies-1)*7+[1:7]) = LL;
        ModelZero.pdoOptions.props.pdoOutputRange = 3*ones(1,length(ModelZero.species));
        if isempty(outputMax)
            ModelZero.pdoOptions.props.pdoOutputRange(obsvSpecies) = ModelZero.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor.size(2);
        else
            ModelZero.pdoOptions.props.pdoOutputRange(obsvSpecies) = outputMax;
        end
        ModelZero.pdoOptions.PDO = ModelZero.generatePDO(ModelZero.pdoOptions,...
            ModelZero.pdoOptions.props.PDOpars,FSPsoln.fsp,false);
end

ModelZero.fittingOptions.timesToFit = zeros(size(ModelZero.dataSet.times),'logical');
for i=1:length(vars.timeSet)
    ModelZero.fittingOptions.timesToFit(ModelZero.dataSet.times==vars.timeSet(i))=true;
end
if strcmp(vars.modelVarsToFit,'all')
    ModelZero.fittingOptions.modelVarsToFit = ones(1,size(ModelZero.parameters,1),'logical');
else
    ModelZero.fittingOptions.modelVarsToFit = vars.modelVarsToFit;
end

end

function [Model,mlikelihood] = runFittingFun(Model,vars,iPDO,pdoPars,...
    muLog10Prior,sigLog10Prior,chainResults,fitOptions,parGuesses)

Model.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior(vars.modelVarsToFit)).^2./(2*sigLog10Prior(vars.modelVarsToFit).^2);

Model.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
Model.initialTime = -1e-6;
Model.fspOptions.usePiecewiseFSP = false;  % Set FSP error tolerance.
Model.fspOptions.initApproxSS = true;  % Set FSP to use SS approximation for IC.

% Solve initial model and compute likelihood
Model.fspOptions.fspTol = 0.000001;  % Set strict FSP error tolerance for current best parameters.
[FSPsoln,Model.fspOptions.bounds] = Model.solve;  % Solve the FSP analysis
Model = associateData(Model,vars,iPDO,pdoPars,FSPsoln);
Model.tSpan = sort(unique([Model.initialTime,5,Model.dataSet.times(Model.fittingOptions.timesToFit)]));
Model.fspOptions.bounds = vars.initialFspBounds;
[FSPsoln,Model.fspOptions.bounds] = Model.solve;  % Solve the FSP analysis
mlikelihood = -Model.computeLikelihood([Model.parameters{Model.fittingOptions.modelVarsToFit,2}]',FSPsoln.stateSpace);
disp(['starting fit found for iPDO = ',num2str(iPDO),'. -logL = ',num2str(mlikelihood)])

% Add best parameter set from MH search.
try
    J = chainResults.mhValue~=0;
    chainResults.mhSamples = chainResults.mhSamples(J,:);
    chainResults.mhValue = chainResults.mhValue(J);
    [~,ir] = max(chainResults.mhValue);
    parGuesses = [exp(chainResults.mhSamples(ir,:));parGuesses];
catch
end

% Run through parGuesses to see if there is a better starting point.
ModTmp = Model;
for jfit=1:size(parGuesses,1)
     x0 = parGuesses(jfit,:);
     ModTmp.parameters(ModTmp.fittingOptions.modelVarsToFit,2) = num2cell(x0);
     ModTmp.tSpan = sort(unique([ModTmp.initialTime,5,ModTmp.dataSet.times(ModTmp.fittingOptions.timesToFit)]));
     ModTmp.fspOptions.bounds = vars.initialFspBounds;
     [FSPsoln,ModTmp.fspOptions.bounds] = ModTmp.solve;  % Solve the FSP analysis
     ModTmp = associateData(ModTmp,vars,iPDO,pdoPars,FSPsoln);
     mlikelihoodGuess = -ModTmp.computeLikelihood([ModTmp.parameters{ModTmp.fittingOptions.modelVarsToFit,2}]',FSPsoln.stateSpace);
     if mlikelihoodGuess<mlikelihood
         mlikelihood = mlikelihoodGuess;
         disp(['improved fit found for iPDO = ',num2str(iPDO),' from ',num2str(jfit),'. -logL = ',num2str(mlikelihood)])
         % Update Model Parameters
         Model = ModTmp;
     end
end

% Here we call the search process with some fitting options.
for ifit=1%:3
    %     % Try fits from other cases.
    %     if ifit==1
    %         modArr = [0,1:iPDO-1,iPDO+1:5];
    %     else
    %         modArr = [0];
    %     end
    %     for jfit = modArr
    %         if jfit==0
    x0 = [Model.parameters{Model.fittingOptions.modelVarsToFit,2}]';
    %         else
    %             x0 = parGuesses(jfit,:);
    %         end

    % Regenerate FSP projection size
    ModTmp = Model;
%     ModTmp.fspOptions.bounds=[];
%     ModTmp.fspOptions.fspTol = 0.00001;  % Set strict FSP error tolerance for current best parameters.
%     ModTmp.parameters(ModTmp.fittingOptions.modelVarsToFit,2) = num2cell(x0);
%     [FSPsoln,ModTmp.fspOptions.bounds] = ModTmp.solve;  % Solve the FSP analysis
%     ModTmp = associateData(ModTmp,vars,iPDO,pdoPars,FSPsoln);

    % Stop expansion during fit
    ModTmp.fspOptions.fspTol = inf;  % Set high FSP error tolerance for faster search.
    [parsMH] = ModTmp.maximizeLikelihood(x0,fitOptions);

    % Check if better fit found during search
    ModTmp.parameters(ModTmp.fittingOptions.modelVarsToFit,2) = num2cell(parsMH);
    ModTmp.fspOptions.bounds = vars.initialFspBounds;
    ModTmp.fspOptions.fspTol = 0.00001;  % Set strict FSP error tolerance for current best parameters.
    [FSPsoln,ModTmp.fspOptions.bounds] = ModTmp.solve;  % Solve the problem
    ModTmp = associateData(ModTmp,vars,iPDO,pdoPars,FSPsoln);
    mlikelihoodGuess = -ModTmp.computeLikelihood(parsMH,FSPsoln.stateSpace);

    if mlikelihoodGuess<mlikelihood
        mlikelihood = mlikelihoodGuess;
        disp(['improved fit found for iPDO = ',num2str(iPDO),'. -logL = ',num2str(mlikelihood)])

        % Update Model Parameters
        Model = ModTmp;
    end
    %     end
end



    % Update Model Parameters
%     Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(pars);

    % Refresh FSP Bounds and PDO for current parameter set.
    %     Model.fspOptions.fspTol = 0.00001;  % Set strict FSP error tolerance for current best parameters.
    %     [FSPsoln,Model.fspOptions.bounds] = Model.solve;  % Solve the FSP analysis
    %     Model = associateData(Model,vars,iPDO,pdoPars,FSPsoln);

    %     Model.makeFitPlot
% end 
Model.fspOptions.fspTol = 0.001;  % Set FSP error tolerance.
Model.dataSet=[]; % Clear data before returning (to save file space).
Model.pdoOptions= struct('unobservedSpecies',[],'PDO',[]); % Options for FIM analyses
end

function [ModelZero,chainResults] = createDefaults(muLog10Prior,vars)
for iPDO = 1:5
    ModelZero{iPDO} = SSIT;    % Create SSIT instance and call it 'Model'.
    switch vars.modelChoice
        case '2stateBurst'
            ModelZero{iPDO}.species = {'x1';'x2';'x3';'x4'}; % Set species names for bursting gene expression model:
            ModelZero{iPDO}.stoichiometry = [-1, 1, 0, 0, 0, 0;...
                1,-1,-1, 1, 0, 0;...
                0, 0, 1,-1, 0, 0;...
                0, 0, 0, 0, 1,-1];
            ModelZero{iPDO}.propensityFunctions = {'k12*x1';'k21*x2';'omega*x2';'1000*x3';'b*1000*x3*ITrypt';'g*x4'};
            ModelZero{iPDO}.initialCondition = [2;0;0;0];
            ModelZero{iPDO}.parameters = ({'k12',0.01;'k21',0.01;'omega',.2;'b',7.12;'g',0.012});
            ModelZero{iPDO}.inputExpressions = {'ITrypt','(t<5)'};
        case '1stateBurst'
            ModelZero{iPDO}.species = {'x2';'x3';'x4'}; % Set species names for bursting gene expression model:
            ModelZero{iPDO}.stoichiometry = [-1, 1, 0, 0;...
                1,-1, 0, 0;...
                0, 0, 1,-1];
            ModelZero{iPDO}.propensityFunctions = {'omega*x2';'1000*x3';'b*1000*x3*ITrypt';'g*x4'};
            ModelZero{iPDO}.initialCondition = [2;0;0];
            ModelZero{iPDO}.parameters = ({'omega',.2;'b',7.12;'g',0.012});
            ModelZero{iPDO}.inputExpressions = {'ITrypt','(t<5)'};
        case '3state'
            ModelZero{iPDO}.species = {'x1';'x2';'x3';'x4'}; % Set species names for bursting gene expression model:
            ModelZero{iPDO}.stoichiometry = [-1, 1, 0, 0, 0, 0;...
                1,-1,-1, 1, 0, 0;...
                0, 0, 1,-1, 0, 0;...
                0, 0, 0, 0, 1,-1];
            ModelZero{iPDO}.propensityFunctions = {'k12*x1';'k21*x2';'omega*x2';'kx*x3';'b*kx*x3*ITrypt';'g*x4'};
            ModelZero{iPDO}.initialCondition = [2;0;0;0];
            ModelZero{iPDO}.parameters = ({'k12',0.01;'k21',0.01;'omega',.2;'b',7.12;'g',0.012;'kx',1000});
            ModelZero{iPDO}.inputExpressions = {'ITrypt','(t<5)'};
        case '3stateWithTS'
%             ModelZero{iPDO}.fspOptions.usePiecewiseFSP=true;
            ModelZero{iPDO}.species = {'x1';'x2';'x3';'x4';'x5'}; % Set species names for bursting gene expression model:
            ModelZero{iPDO}.stoichiometry = [-1, 1, 0, 0, 0, 0, 0;...
                1,-1,-1, 1, 0, 0, 0;...
                0, 0, 1,-1, 0, 0, 0;...
                0, 0, 0, 0, 0, 1, -1;...
                0, 0, 0, 0, 1,-1, 0];
            ModelZero{iPDO}.propensityFunctions = {'k12*x1';'k21*x2';'omega*x2';'kx*x3';'b*kx*x3*ITrypt';'kesc*x5';'g*x4'};
            ModelZero{iPDO}.initialCondition = [2;0;0;0;0];
            ModelZero{iPDO}.parameters = ({'k12',0.01;'k21',0.01;'omega',.2;'b',7.12;'g',0.012;'kx',1000;'kesc',0.2});
            ModelZero{iPDO}.inputExpressions = {'ITrypt','(t<0)'};
        case '2statePoisson'
%             ModelZero{iPDO}.fspOptions.usePiecewiseFSP=true;
            ModelZero{iPDO}.species = {'x1';'x2';'x4'}; % Set species names for bursting gene expression model:
            ModelZero{iPDO}.stoichiometry = [-1, 1, 0, 0;...
                1,-1, 0, 0;...
                0, 0, 1,-1];
            ModelZero{iPDO}.propensityFunctions = {'k12*x1';'k21*x2';'omega*x2*ITrypt';'g*x4'};
            ModelZero{iPDO}.initialCondition = [2;0;0];
            ModelZero{iPDO}.parameters = ({'k12',0.01;'k21',0.01;'omega',.2;'g',0.012});
            ModelZero{iPDO}.inputExpressions = {'ITrypt','(t<5)'};
    end

    ModelZero{iPDO}.parameters(:,2) = num2cell(10.^muLog10Prior);
    chainResults{iPDO}=[];
end
end

function [ModelZero,chainResults] = createModelDusp1(muLog10Prior)
for iPDO = 5:-1:1
ModelZero = SSIT;
ModelZero.species = {'x1';'x2'};
ModelZero.initialCondition = [0;0];
ModelZero.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
ModelZero.stoichiometry = [1,-1,0,0;0,0,1,-1];
ModelZero.inputExpressions = {'IGR','a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'};
ModelZero.parameters = ({'koff',0.14;'kon',0.14;'kr',25;'gr',0.01;...
    'a0',0.006;'a1',0.4;'r1',0.04;'r2',0.1});
ModelZero.initialTime = -1e-6;  % t0<<0 to simulate steady state at t=0
ModelZero.fspOptions.initApproxSS = true;  % Set FSP to use SS approximation for IC.
%%
ModelZero{iPDO}.parameters(:,2) = num2cell(10.^muLog10Prior);
chainResults{iPDO}=[];

end
end

function makePatch(xvals,logmnIdeal,logstdIdeal,col)

X = [xvals,xvals(end:-1:1)];
Y = [logmnIdeal-logstdIdeal,logmnIdeal(end:-1:1)+logstdIdeal(end:-1:1)];
patch(X,Y,col)

end

