function results = FittingFunctionsCoLocalized(iStep,fileName,vars)
arguments
    iStep
    fileName
    vars = [];
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
    'pdoTimes',[0];...
    'FIMnCells',[0,0,0];...
    'dataFile','Huy_paper_data/SpotClassification_ts0_550_ts1_400dTS_2.csv'};

for i=1:size(defaults,1)
    if ~isfield(vars,defaults{i,1})
        vars.(defaults{i,1}) = defaults{i,2};
    end
end
binSize = 20;
chainResults=cell(1,5);
fitOptions = optimset('Display',vars.display,'MaxIter',vars.iter);
muLog10Prior = [-4,-4,log10(0.2),log10(7.12),log10(0.012)]';
sigLog10Prior = [1,1,.5,.5,.5]'*vars.priorScale;

%%
switch iStep
    case 1
        try
            load(fileName,'lTrue2FISH','lTrue2MCP')
        catch
            lTrue2FISH = [1,0.5,1.0];
            lTrue2MCP = [1,0.5,1.0];
        end

        DATA = importdata(vars.dataFile);        
        % 0 = MCP-GFP
        % 1 = smiFISH
        J = find(strcmp(DATA.colheaders,'number_spots_type_1'));
        K = find(strcmp(DATA.colheaders,'number_spots_type_0'));
        I = find(strcmp(DATA.colheaders,'total_spots'));
        jT = find(strcmp(DATA.colheaders,'time'));

        jCells = DATA.data(:,jT)==vars.pdoTimes(1);
        for i=2:length(vars.pdoTimes)
            jCells(DATA.data(:,jT)==vars.pdoTimes(i))=true;
        end

        DistortedFISH = DATA.data(jCells,J);
        DistortedMCP = DATA.data(jCells,K);
        True = DATA.data(jCells,I);

        if vars.doFit==1
            OBJ = @(l)-findError(l,True,DistortedFISH);
            for i=1:1
                lTrue2FISH = fminsearch(OBJ,lTrue2FISH,fitOptions);
            end

            OBJ2 = @(l)-findError(l,True,DistortedMCP);
            for i=1:1
                lTrue2MCP = fminsearch(OBJ2,lTrue2MCP,fitOptions);
            end

            try
                save(fileName,'lTrue2FISH','lTrue2MCP','-append')
            catch
                save(fileName,'lTrue2FISH','lTrue2MCP')
            end
            lTrue2FISH
            lTrue2MCP
        else

            [~,cTrue2FISH] = findError(lTrue2FISH,True,DistortedFISH);
            [~,cTrue2MCP] = findError(lTrue2MCP,True,DistortedMCP);

            close all

            for icase = 1:3
                switch icase
                    case 1
                        jCells = DATA.data(:,jT)==0;
                    case 2
                        jCells = DATA.data(:,jT)==18;
                    case 3
                        jCells = DATA.data(:,jT)==300;
                end

                DistortedFISH = DATA.data(jCells,J);
                DistortedMCP = DATA.data(jCells,K);
                True = DATA.data(jCells,I);

                %%
                figure(2)
                switch icase
                    case 1
                        Z = max(-200,log10(cTrue2FISH));
                        contourf([0:size(cTrue2FISH,2)-1],[0:size(cTrue2FISH,1)-1],Z)
                        colorbar
                        hold on
                        scatter(True,DistortedFISH,100,'sk','filled')
                    case 2
                        scatter(True,DistortedFISH,100,'om','filled')
                    case 3
                        scatter(True,DistortedFISH,100,'^c','filled')
                end

                xlabel('Number observed (total)')
                ylabel('Number observed (smiFISH)')
                set(gca,'fontsize',16,'xlim',[0,700],'ylim',[0,700])

                %%
                figure(22)
                switch icase
                    case 1
                        Z = max(-200,log10(cTrue2MCP));
                        contourf([0:size(cTrue2MCP,2)-1],[0:size(cTrue2MCP,1)-1],Z)
                        colorbar
                        hold on
                        scatter(True,DistortedMCP,100,'sk','filled')
                    case 2
                        scatter(True,DistortedMCP,100,'om','filled')
                    case 3
                        scatter(True,DistortedMCP,100,'^c','filled')
                end

                xlabel('Number observed (total)')
                ylabel('Number observed (MCP-GFP)')
                set(gca,'fontsize',16,'xlim',[-1,700],'ylim',[-1,700])

                col=[0.2,0.2,0.6;...
                    0.6,0.2,0.2;...
                    0.2,0.6,0.6;...
                    0.4,0.4,0.6;...
                    0.6,0.4,0.4;...
                    0.4,0.6,0.4;...
                    0.6,0.6,0.9;...
                    0.9,0.6,0.6;...
                    0.6,0.9,0.6];

                %%
                figure(10);
                Nmax = min(1000,size(cTrue2FISH,2));
                edges = -1:Nmax;
                h1 = histogram(DistortedFISH,edges,'normalization','cdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(1,:)); hold on
                h2 = histogram(True,edges,'normalization','cdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(2,:));

                figure(5+100*icase)
                cdf1 = h2.Values;
                cdf2 = h1.Values;
                stairs([-1:Nmax-1],h2.Values,'LineWidth',4,'Color',col(2+3*(icase-1),:)); hold on
                stairs([-1:Nmax-1],h1.Values,'LineWidth',4,'Color',col(1+3*(icase-1),:)); hold on
                close(10)

                set(gca,'FontSize',16,'xlim',[-1,800],'ylim',[0,1.02])
                xlabel('Number of mRNA')
                ylabel('Cumulative Probability')

                figure(6+100*icase);
                binSize = 50;
                edges = 0:binSize:Nmax;
                histogram(DistortedFISH,edges,'normalization','pdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(1,:)); hold on
                histogram(True,edges,'normalization','pdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(2,:));
                set(gca,'FontSize',16,'xlim',[-1,800])
                xlabel('Number of mRNA')
                ylabel('Probability Mass')

                figure(10); clf
                edges = 0:Nmax;
                h1 = histogram(DistortedFISH,edges,'normalization','pdf');
                h1 = h1.Values;
                h2 = histogram(True,edges,'normalization','pdf');
                h2 = h2.Values;
                close(10);

                pDist = cTrue2FISH(:,1:Nmax);
                distPred = pDist*h2';

                pDistBin=[];
                for i = 1:floor(length(distPred)/binSize)
                    pDistBin(i) = sum(distPred((i-1)*binSize+1:i*binSize))/binSize;
                end
                figure(5+100*icase)
                stairs([-1:length(distPred)],[0;cumsum(distPred);1],'--','LineWidth',4,'Color',col(3,:));

                KS = max(abs([0;cumsum(distPred);ones(length(cdf2)-length(distPred)-1,1)] - cdf2'))
                figure(6+100*icase)
                stairs([0:binSize:length(distPred)-binSize],(pDistBin),'LineWidth',4,'Color',col(3,:));
                legend('smFISH' ,'total' ,'total + PDO')

                %%
                figure(20);
                Nmax = min(1000,size(cTrue2MCP,2));
                edges = -1:Nmax;
                h1 = histogram(DistortedMCP,edges,'normalization','cdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(1,:)); hold on
                h2 = histogram(True,edges,'normalization','cdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(2,:));
                figure(25+100*icase)
                cdf1 = h2.Values;
                cdf2 = h1.Values;
                stairs([-1:Nmax-1],h2.Values,'LineWidth',4,'Color',col(2,:)); hold on
                stairs([-1:Nmax-1],h1.Values,'LineWidth',4,'Color',col(1,:)); hold on
                close(20)

                set(gca,'FontSize',16,'xlim',[-1,800],'ylim',[0,1.02])
                xlabel('Number of mRNA')
                ylabel('Cumulative Probability')

                figure(26+100*icase);
                edges = 0:binSize:Nmax;
                histogram(DistortedMCP,edges,'normalization','pdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(1,:)); hold on
                histogram(True,edges,'normalization','pdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col(2,:));
                set(gca,'FontSize',16,'xlim',[0,800])
                xlabel('Number of mRNA')
                ylabel('Probability Mass')

                figure(20); clf
                edges = 0:Nmax;
                h1 = histogram(DistortedMCP,edges,'normalization','pdf');
                h1 = h1.Values;
                h2 = histogram(True,edges,'normalization','pdf');
                h2 = h2.Values;
                close(20);

                pDist = cTrue2MCP(:,1:Nmax);
                distPred = pDist*h2';

                pDistBin=[];
                for i = 1:floor(length(distPred)/binSize)
                    pDistBin(i) = sum(distPred((i-1)*binSize+1:i*binSize))/binSize;
                end
                figure(25+100*icase)
                stairs([-1:length(distPred)],[0;cumsum(distPred);1],'--','LineWidth',4,'Color',col(3,:));

                KS = max(abs([0;cumsum(distPred);ones(length(cdf2)-length(distPred)-1,1)] - cdf2'))
                figure(26+100*icase)
                stairs([0:binSize:length(distPred)-binSize],(pDistBin),'LineWidth',4,'Color',col(3,:));
                legend('MCP-GFP' ,'total' ,'total + PDO')

                figure(5+100*icase)
                legend('total' ,'smiFISH' ,'total + PDO')

                figure(2+100*icase)
                legend('PDO contours','data')

                figure(25+100*icase)
                legend('total' ,'MCP-GFP' ,'total + PDO')

                figure(22+100*icase)
                legend('PDO contours','data')

            end
        end

    case 32 % Fit model to data
        % For this fit, we will assume that the smFISH data is "true" and we will
        % fit once using the smFISH data and once with the MCP-GFP data.

        load(fileName,'ModelZero','lTrue2FISH','lTrue2MCP','chainResults')
        try
            if ~exist('ModelZero','var')
                error('missing ModelZero')
            end
        catch
            disp('Could not find previous model.  Starting from scratch.')
            [ModelZero,chainResults] = createDefaults(muLog10Prior);
        end

        if isempty(vars)||vars.doFit==1
            parGuessesLocal =[];
            for iPDO = 1:5
                parGuessesLocal = [parGuessesLocal;[ModelZero{iPDO}.parameters{vars.modelVarsToFit,2}]];
            end

            parfor iPDO = 1:5
                ModelZero{iPDO}.fspOptions.usePiecewiseFSP=false;
                ModelZero{iPDO}.fspOptions.initApproxSS=true;
                ModelZero{iPDO}.inputExpressions = {'ITrypt','(t<5)'};
                ModelZero{iPDO} = runFittingFun(ModelZero{iPDO},vars,iPDO,...
                    lTrue2FISH,lTrue2MCP,muLog10Prior,sigLog10Prior,chainResults{iPDO},...
                    fitOptions,parGuessesLocal);
            end
            save(fileName,'ModelZero','-append')

        else

            parsTrue = ModelZero{1}.parameters;
            parsDistortedMCP = ModelZero{2}.parameters;
            parsDistortedFISH = ModelZero{3}.parameters;
            parsCorrectedMCP = ModelZero{4}.parameters;
            parsCorrectedFISH = ModelZero{5}.parameters;
            T = table([parsTrue{:,2}]',[parsDistortedMCP{:,2}]',[parsDistortedFISH{:,2}]',[parsCorrectedMCP{:,2}]',[parsCorrectedFISH{:,2}]')

            diffTrueDistortedMCP = exp(abs(log([parsDistortedMCP{:,2}]./[parsTrue{:,2}])))';
            diffTrueDistortedFISH = exp(abs(log([parsDistortedFISH{:,2}]./[parsTrue{:,2}])))';
            diffTrueCorrectedMCP = exp(abs(log([parsCorrectedMCP{:,2}]./[parsTrue{:,2}])))';
            diffTrueCorrectedFISH = exp(abs(log([parsCorrectedFISH{:,2}]./[parsTrue{:,2}])))';

            meanFoldDiffDistMCP = mean(diffTrueDistortedMCP(ModelZero{1}.fittingOptions.modelVarsToFit))
            meanFoldDiffCorectedMCP = mean(diffTrueCorrectedMCP(ModelZero{1}.fittingOptions.modelVarsToFit))
            meanFoldDiffDistFISH = mean(diffTrueDistortedFISH(ModelZero{1}.fittingOptions.modelVarsToFit))
            meanFoldDiffCorectedFISH = mean(diffTrueCorrectedFISH(ModelZero{1}.fittingOptions.modelVarsToFit))

            for iPDO = 1:5
                ModelZero{iPDO}.inputExpressions = {'ITrypt','(t<5)'};
                ModelZero{iPDO}.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
                ModelZero{iPDO}.fspOptions.fspTol = 0.0001;  % Set FSP error tolerance.
                ModelZero{iPDO}.tSpan = [-1e-6,0,5,18,300];

                [FSPsoln,ModelZero{iPDO}.fspOptions.bounds] = ModelZero{iPDO}.solve;  % Solve the FSP analysis
                ModelZero{iPDO} = associateData(ModelZero{iPDO},vars,iPDO,lTrue2FISH,lTrue2MCP,FSPsoln);
                ModelZero{iPDO}.tSpan = sort(unique([ModelZero{iPDO}.initialTime,5,ModelZero{iPDO}.dataSet.times]));

                [~,~,fitSolution{iPDO}] = ModelZero{iPDO}.computeLikelihood;

                results.fitLikelihoods(iPDO,:) = (fitSolution{iPDO}.DataLoadingAndFittingTabOutputs.V_LogLk-...
                    fitSolution{iPDO}.DataLoadingAndFittingTabOutputs.perfectMod);
            end

            ModelTMP1 = ModelZero{1};
            ModelTMP1.parameters = ModelZero{4}.parameters;
            [~,~,fitSolution{6}] = ModelTMP1.computeLikelihood;
            results.fitLikelihoods(6,:) = (fitSolution{6}.DataLoadingAndFittingTabOutputs.V_LogLk-...
                fitSolution{6}.DataLoadingAndFittingTabOutputs.perfectMod);

            ModelTMP = ModelZero{1};
            ModelTMP.parameters = ModelZero{5}.parameters;
            [~,~,fitSolution{7}] = ModelTMP.computeLikelihood;
            results.fitLikelihoods(7,:) = (fitSolution{7}.DataLoadingAndFittingTabOutputs.V_LogLk-...
                fitSolution{7}.DataLoadingAndFittingTabOutputs.perfectMod);

            results.fitSolution=fitSolution;

            ModelZero{1}.makeFitPlot(fitSolution{1},binSize,[1,2,3,4])  % Fit total
            ModelZero{1}.makeFitPlot(fitSolution{1},binSize,[10,20,3,4]) % Fit total
            ModelZero{2}.makeFitPlot(fitSolution{2},binSize,[1,2,3,4]) % Fit MCP no correction
            ModelZero{3}.makeFitPlot(fitSolution{3},binSize,[10,20,3,4]) % Fit smiFISH no correction
            ModelZero{4}.makeFitPlot(fitSolution{4},binSize,[1,2,3,4]) % Fit MCP with correction
            ModelZero{5}.makeFitPlot(fitSolution{5},binSize,[10,20,3,4]) % Fit smiFISH with correction

            ModelTMP1.makeFitPlot(fitSolution{6},binSize,[1,2,3,4]);
            ModelTMP.makeFitPlot(fitSolution{7},binSize,[10,20,3,4]);

            for ifig=[1,2,10,20]
                figure(ifig)
                for isub = 1:3
                    subplot(1,3,isub)
                    h = gca;
                    % Total mRNA - black
                    h.Children(end).Color = [0,0,0];
                    h.Children(end-1).Color = [0.5,0.5,0.5];
                    %             h.Children(end-2).Color = [0,0,0];
                    h.Children(2).Color = [0,0,0];
                    % MCP data
                    h.Children(end-2).Color = 'b';
                    h.Children(end-4).Color = 'b';

                    % MCP un-corrected fit
                    h.Children(5).LineWidth = 4;
                    h.Children(5).Color = 'm';

                    % MCP corrected fit
                    h.Children(3).LineWidth = 3;
                    h.Children(3).Color = 'c';

                    % prediction of total
                    h.Children(1).Color = 'c';
                    h.Children(1).LineStyle = '--';

                end
            end
        end

    case 3 % Compute FIM for best parameter fit
        load(fileName,'ModelZero','lTrue2FISH','lTrue2MCP')
        for iPDO = 1:5
            ModelZero{iPDO}.fspOptions.fspTol = 1e-6;  % Set FSP error tolerance.
            ModelZero{iPDO}.solutionScheme = 'FSP'; % Set solutions scheme to FSP
            [FSPsoln,ModelZero{iPDO}.fspOptions.bounds] = ModelZero{iPDO}.solve;  % Solve the problem
            ModelZero{iPDO} = associateData(ModelZero{iPDO},vars,iPDO,lTrue2FISH,lTrue2MCP,FSPsoln);

            ModelZero{iPDO}.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
            ModelZero{iPDO}.sensOptions.solutionMethod = 'finiteDifference';

            ModelZero{iPDO}.tSpan = unique([[-1e-6,5],...
                ModelZero{iPDO}.dataSet.times(ModelZero{iPDO}.fittingOptions.timesToFit)]);

            [sensSoln,ModelZero{iPDO}.fspOptions.bounds] = ModelZero{iPDO}.solve;  % Solve the sensitivity problem

            try
                fimResults{iPDO} = ModelZero{iPDO}.computeFIM(sensSoln.sens);
            catch
                if iPDO>=4
                    ModelZero{iPDO}.fspOptions.fspTol = 1e-7;  % Set FSP error tolerance.
                    ModelZero{iPDO}.solutionScheme = 'FSP'; % Set solutions scheme to FSP
                    [FSPsoln,ModelZero{iPDO}.fspOptions.bounds] = ModelZero{iPDO}.solve;  % Solve the problem
                    ModelZero{iPDO}.pdoOptions.PDO = ModelZero{iPDO}.generatePDO(ModelZero{iPDO}.pdoOptions,...
                        ModelZero{iPDO}.pdoOptions.props.PDOpars,FSPsoln.fsp,false);
                end
                fimResults{iPDO} = ModelZero{iPDO}.computeFIM(sensSoln.sens);
            end

            cellCounts = zeros(1,length(ModelZero{iPDO}.tSpan));
            for it = 1:length(ModelZero{iPDO}.dataSet.times)
                if ModelZero{iPDO}.fittingOptions.timesToFit(it)
                    [~,jm] = min(abs(ModelZero{iPDO}.tSpan-ModelZero{iPDO}.dataSet.times(it)));
                    cellCounts(jm) = ModelZero{iPDO}.dataSet.nCells(it);
                end
            end

            FIMZero{iPDO} = ModelZero{iPDO}.evaluateExperiment(fimResults{iPDO},cellCounts);

            FIMZeroLog{iPDO} = diag([ModelZero{iPDO}.parameters{ModelZero{iPDO}.fittingOptions.modelVarsToFit,2}])*...
                FIMZero{iPDO}(ModelZero{iPDO}.fittingOptions.modelVarsToFit,ModelZero{iPDO}.fittingOptions.modelVarsToFit)*...
                diag([ModelZero{iPDO}.parameters{ModelZero{iPDO}.fittingOptions.modelVarsToFit,2}]);
            covLog{iPDO} = FIMZeroLog{iPDO}^-1;
        end
        save(fileName,'FIMZero','FIMZeroLog','fimResults','covLog','-append')

    case 4 % Run Met Hast to quantify uncertainty on initial parameters
        if isempty(vars)||vars.doFit==1
            load(fileName,'ModelZero','FIMZeroLog')
            for iPDO = 1:5
                indsPars = ModelZero{iPDO}.fittingOptions.modelVarsToFit;
                ModelZero{iPDO}.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior(indsPars)).^2./(2*sigLog10Prior(indsPars).^2);
            end
            for iPDO = 1:5
                % Add a small perturbation to log-FIM to make sure that it will
                % be invertible.
                if cond(FIMZeroLog{iPDO})>1e7
                    covLog = (FIMZeroLog{iPDO} + 1e-2*eye(length(FIMZeroLog{iPDO})))^-1;
                else
                    covLog = (FIMZeroLog{iPDO})^-1;
                end

                MHOptions = struct('numberOfSamples',vars.nMH,'burnIn',0,'thin',3,'saveFile',['TMPmh_',num2str(iPDO),'_',fileName]);
                proposalWidthScale = vars.mhScaling;
                MHOptions.proposalDistribution  = @(x)mvnrnd(x,proposalWidthScale*(covLog+covLog')/2);

                ModelZero{iPDO}.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
                ModelZero{iPDO}.fspOptions.fspTol = 0.0001;  % Set FSP error tolerance.
                [~,ModelZero{iPDO}.fspOptions.bounds] = ModelZero{iPDO}.solve;  % Solve the FSP analysis

                ModelZero{iPDO}.fspOptions.fspTol = inf;  % Set FSP error tolerance.

                parGuess = [ModelZero{iPDO}.parameters{ModelZero{iPDO}.fittingOptions.modelVarsToFit,2}];
                [pars,~,chainResults{iPDO}] = ModelZero{iPDO}.maximizeLikelihood(parGuess',MHOptions,'MetropolisHastings');

                ModelZero{iPDO}.fspOptions.fspTol = 0.001;  % Set FSP error tolerance.

                ModelZero{iPDO}.parameters(ModelZero{iPDO}.fittingOptions.modelVarsToFit,2) =...
                    num2cell(pars(1:length(ModelZero{iPDO}.fittingOptions.modelVarsToFit)));

                %                 parGuess = [ModelZero{iPDO}.parameters{ModelZero{iPDO}.fittingOptions.modelVarsToFit,2}]';
                %                 [~,~,fitSolutions] = ModelZero{iPDO}.computeLikelihood(parGuess);
            end
            save(fileName,'ModelZero','chainResults','-append')
            for i=1:5
                disp('MH Completed successfully. Deleting temporary files.')
                delete(['TMPmh_',num2str(i),'_',fileName,'.mat'])
            end

        elseif vars.doFit==0
            load(fileName,'chainResults','covLog','ModelZero')
            for iPDO = 1:5
                log10MLEPars = log10([ModelZero{iPDO}.parameters{ModelZero{iPDO}.fittingOptions.modelVarsToFit,2}]);
                smplDone = chainResults{iPDO}.mhSamples(:,:);
                valDone = chainResults{iPDO}.mhValue;
                smplDone = smplDone(valDone~=0,:);
                valDone = valDone(valDone~=0);

                %                 N = floor(length(valDone)/2);
                %                 smplDone = smplDone(N+1:end,:);
                %                 valDone = valDone(N+1:end);

                [valDoneSorted,J] = sort(valDone);
                smplDone = smplDone(J,:);

                covPRIOR = diag(sigLog10Prior(ModelZero{iPDO}.fittingOptions.modelVarsToFit));
                muPRIOR = muLog10Prior(ModelZero{iPDO}.fittingOptions.modelVarsToFit);
                covFIM = covLog{iPDO}/(log(10)^2);
                covMLE = cov(smplDone/log(10));
                stdMLE = sqrt(diag(covMLE))';
                stdFIM = sqrt(diag(covFIM))';
                stdsFIM(iPDO,:) = stdFIM;
                stdsMLE(iPDO,:) = stdMLE;
                detsCovFIM(iPDO) = det(covFIM);
                detsCovMLE(iPDO) = det(covMLE);

                figure(123);
                subplot(1,5,iPDO);

                if nargout>=1
                    results.covMLE{iPDO} = covMLE;
                    results.detCovMLE = detsCovMLE;
                    results.muMLE(iPDO,:) = mean(exp(smplDone));
                    results.MLE(iPDO,:) = log10MLEPars;
                    results.sigMLE(iPDO,:) = std(exp(smplDone));
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

                            ssit.parest.ellipse(log10MLEPars([j,i]),icdf('chi2',0.9,2)*covLog{iPDO}([j,i],[j,i])/(log(10)^2),'k-','linewidth',3)
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
            figure
            bar([stdsFIM(1,:);stdsMLE(1,:);stdsFIM(4,:);stdsMLE(4,:);stdsFIM(5,:);stdsMLE(5,:)]')
            figure
            bar([detsCovFIM([1,4,5]);detsCovMLE([1,4,5])]')
        end
    case 5 % Sample from MH Results and Compute FIM at all time points
        clear ModelOne sensSolnOne fimResultsOne
        load(fileName,'ModelZero','chainResults','lTrue2FISH','lTrue2MCP')
        nSamples = vars.nFIMsamples;

        parfor iPDO = 1:5
            smplDone = chainResults{iPDO}.mhSamples(:,:);
            ismpls = floor(linspace(1,size(smplDone,1),nSamples));
            parsets{iPDO} = repmat([ModelZero{iPDO}.parameters{:,2}],nSamples,1);
            parsets{iPDO}(:,vars.modelVarsToFit) = exp(smplDone(ismpls,:));
            parsets{iPDO}(1,:) = [ModelZero{iPDO}.parameters{:,2}];

            for i=1:nSamples
                [i,iPDO]

                ModelOne = ModelZero{iPDO};
                ModelOne.tSpan = sort(unique([-1e-6,0:6:1200,5]));
                ModelOne.parameters(:,2) = num2cell(parsets{iPDO}(i,:));

                ModelOne.solutionScheme = 'FSP'; % Set solutions scheme to FSP Sensitivity
                ModelOne.fspOptions.fspTol = 1e-6;
                [FSPsoln,ModelOne.fspOptions.bounds] = ModelOne.solve;  % Solve the FSP analysis
                ModelOne = associateData(ModelOne,vars,iPDO,lTrue2FISH,lTrue2MCP,FSPsoln);

                ModelOne.fspOptions.fspTol = 1e-5;
                ModelOne.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
                ModelOne.sensOptions.solutionMethod = 'finiteDifference';
                sensSolnOne = ModelOne.solve;  % Solve the sensitivity problem

                try
                    fimResultsOne{i,iPDO} = ModelOne.computeFIM(sensSolnOne.sens);
                catch
                    if iPDO>=4
                        ModelOne.fspOptions.fspTol = 1e-7;  % Set FSP error tolerance.
                        ModelOne.solutionScheme = 'FSP'; % Set solutions scheme to FSP
                        [FSPsoln,ModelOne.fspOptions.bounds] = ModelOne.solve;  % Solve the problem
                        ModelOne.pdoOptions.PDO = ModelOne.generatePDO(ModelOne.pdoOptions,...
                            ModelOne.pdoOptions.props.PDOpars,FSPsoln.fsp,false);
                    end
                    fimResultsOne{i,iPDO} = ModelOne.computeFIM(sensSolnOne.sens);
                end
            end
        end
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
            for iPDO = 5:-1:1
                for it=size(fimResultsOne{1,1},1):-1:1
                    fimResultsRed{jpar,iPDO}{it,1} = ...
                        diag(parsets{iPDO}(jpar,vars.modelVarsToFit))*...
                        fimResultsOne{jpar,iPDO}{it}(vars.modelVarsToFit,vars.modelVarsToFit)*...
                        diag(parsets{iPDO}(jpar,vars.modelVarsToFit));
                end
            end
        end

        %         NCells = [155 106 66];

        for jpar = 1:size(fimResultsRed,1)
            for iPDO = 1:5
                for it = 2:size(fimResultsRed{1,1},1)
                    nCellsDesign = zeros(1,size(fimResultsRed{1,1},1));
                    nCellsDesign(2) = vars.FIMnCells(1);
                    nCellsDesign(6) = vars.FIMnCells(2);
                    nCellsDesign(53) = vars.FIMnCells(3);
                    nCellsDesign(it) = nCellsDesign(it)+100;
                    [~,covAll,fimMetrics] = ModelZero{1}.evaluateExperiment(fimResultsRed{jpar,iPDO},nCellsDesign);
                    DetDesign(jpar,iPDO,it) = fimMetrics.det;
                    DetCovDesign(jpar,iPDO,it) = det(covAll);
                    EDesign(jpar,iPDO,it) = fimMetrics.minEigVal;
                    TDesign(jpar,iPDO,it) = fimMetrics.trace;
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
        logFIMCritPDOsmiFISH = log10(squeeze(design(:,5,2:end)));
        logFIMCritPDOMCP = log10(squeeze(design(:,4,2:end)));
        logmnPDOFISH = mean(logFIMCritPDOsmiFISH);
        logstdPDOFISH = std(logFIMCritPDOsmiFISH);
        logmnPDOMCP= mean(logFIMCritPDOMCP);
        logstdPDOMCP = std(logFIMCritPDOMCP);

        if length(logFIMCritIDeal)==202
            xvals = sort(unique([0:6:1200,5]));
        else
            xvals = sort(unique([-1e-6,0:6:1200,5]));
        end

        figure

        if ~isfield(vars,'FIMShowErrorBars')||vars.FIMShowErrorBars
            errorbar(xvals,logmnIdeal,logstdIdeal,...
                'color',[0.8,0.4,0.8],'LineWidth',3); hold on
            errorbar(xvals,logmnPDOMCP,logstdPDOMCP,...
                'color',[0.4,0.8,0.4],'LineWidth',3);
            errorbar(xvals,logmnPDOFISH,logstdPDOFISH,...
                'color',[0.8,0.8,0.4],'LineWidth',3);
            plot(xvals,logmnIdeal,...
                'color','k','LineWidth',3); hold on
            plot(xvals,logmnPDOMCP,...
                'color','k','LineWidth',3); hold on
            plot(xvals,logmnPDOFISH,...
                'color','k','LineWidth',3); hold on
        else
            plot(xvals,logmnIdeal,...
                'color',[0.8,0.4,0.8],'LineWidth',3); hold on
            plot(xvals,logmnPDOMCP,...
                'color',[0.4,0.8,0.4],'LineWidth',3);
            plot(xvals,logmnPDOFISH,...
                'color',[0.8,0.8,0.4],'LineWidth',3);
        end


        ylim = get(gca,'ylim');
        tt = xvals;
        if minormax==1
            [~,j1] = min(logmnIdeal);
            [~,j2] = min(logmnPDOFISH);
            [~,j3] = min(logmnPDOMCP);
        else
            [~,j1] = max(logmnIdeal);
            [~,j2] = max(logmnPDOFISH);
            [~,j3] = max(logmnPDOMCP);
        end
        plot(tt(j1)*[1,1],ylim+[-10,10],'color',[0.8,0.4,0.8],'LineWidth',3)
        plot(tt(j3)*[1,1],ylim+[-10,10],'color',[0.4,0.8,0.4],'LineWidth',3)
        plot(tt(j2)*[1,1],ylim+[-10,10],'color',[0.8,0.8,0.4],'LineWidth',3)
        plot(tt(j1)*[1,1],ylim+[-10,10],'k--','LineWidth',3)
        plot(tt(j3)*[1,1],ylim+[-10,10],'k--','LineWidth',3)
        plot(tt(j2)*[1,1],ylim+[-10,10],'k--','LineWidth',3)

        set(gca,'FontSize',16); xlabel('Time (min)'); ylabel(titltxt)
        legend('Ideal','MCP-GFP+PDO','smiFISH+PDO')

end
end

%%  Functions
function [logL,P] = findError(lambda,True,Distorted)
arguments
    lambda
    True
    Distorted
end
% Computes likelihood of observed data given the model of affine poisson
% extra spot counting and probability of measurmeent failure.
Nmax = max(max(max(True,Distorted)));

Np = max(1,lambda(2)+lambda(3)*Nmax);
Np = max(Nmax,ceil(Np+5*sqrt(Np)));
P = zeros(Np+1,Nmax+1);

for xi = 0:Nmax
    % Affine Poisson Gain followed by binomial loss
    P(1:Np+1,xi+1) = pdf('poiss',[0:Np]',max(lambda(1),lambda(2)+lambda(3)*xi));
end


logP = max(log(P),-100);
logL = 0;
for i = 1:length(True)
    logL = logL + logP(Distorted(i)+1,True(i)+1);
end
logL = logL-1e4*(lambda(1)<0);
end

function ModelZero = associateData(ModelZero,vars,iPDO,lTrue2FISH,lTrue2MCP,FSPsoln)

switch iPDO
    case 1 % total spots, no distortion model
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','total_spots'});
        ModelZero.pdoOptions.PDO = [];
    case 2 % MCP-GFP data, no distortion model
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','number_spots_type_0'});
        ModelZero.pdoOptions.PDO = [];
    case 3 % smiFISH data, no distortion model
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','number_spots_type_1'});
        ModelZero.pdoOptions.PDO = [];
    case 4 % MCP-GFP data, Binomial distortion model
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','number_spots_type_0'});
        LL = lTrue2MCP;
        ModelZero.pdoOptions.type = 'AffinePoiss';
        ModelZero.pdoOptions.props.PDOpars = [0,0,0,0,0,0,0,0,0,LL];
        ModelZero.pdoOptions.PDO = ModelZero.generatePDO(ModelZero.pdoOptions,...
            ModelZero.pdoOptions.props.PDOpars,FSPsoln.fsp,false);
    case 5 % smiFISH data, Binomial distortion model
        ModelZero = ModelZero.loadData(vars.dataFile,{'x4','number_spots_type_1'});
        LL = lTrue2FISH;
        ModelZero.pdoOptions.type = 'AffinePoiss';
        ModelZero.pdoOptions.props.PDOpars = [0,0,0,0,0,0,0,0,0,LL];
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
function Model = runFittingFun(Model,vars,iPDO,lTrue2FISH,lTrue2MCP,...
    muLog10Prior,sigLog10Prior,chainResults,fitOptions,parGuesses)

Model.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior(vars.modelVarsToFit)).^2./(2*sigLog10Prior(vars.modelVarsToFit).^2);

Model.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
Model.initialTime = -1e-6;
Model.fspOptions.usePiecewiseFSP = false;  % Set FSP error tolerance.
Model.fspOptions.initApproxSS = true;  % Set FSP to use SS approximation for IC.
Model.tSpan = sort(unique([Model.initialTime,5,Model.dataSet.times(Model.fittingOptions.timesToFit)]));
Model.fittingOptions.pdoVarsToFit = [];

% Here we use the current parameters as our initial guess:
% Here we call the search process with some fitting options.
for ifit=1:3
    % Refresh FSP Bounds and PDO for current parameter set.
    Model.fspOptions.fspTol = 0.00001;  % Set strict FSP error tolerance for current best parameters.
    [FSPsoln,Model.fspOptions.bounds] = Model.solve;  % Solve the FSP analysis
    Model = associateData(Model,vars,iPDO,lTrue2FISH,lTrue2MCP,FSPsoln);
    Model.fspOptions.fspTol = inf;  % Set high FSP error tolerance during search
    
    x0 = [Model.parameters{Model.fittingOptions.modelVarsToFit,2}]';
    [pars,mlikelihood] = Model.maximizeLikelihood(x0,fitOptions);
    % Try fits from other cases.
    if ifit==1
        for jfit = 1:size(parGuesses,1)
            if jfit~=iPDO
                x0 = parGuesses(jfit,:);
                [parsMH,mlikelihoodGuess] = Model.maximizeLikelihood(x0,fitOptions);
                if mlikelihoodGuess<mlikelihood
                    disp(['improved fit found for iPDO = ',num2str(iPDO)])
                    mlikelihood = mlikelihoodGuess
                    pars = parsMH
                end
            end
        end
    end

    try
        J = chainResults.mhValue~=0;
        chainResults.mhSamples = chainResults.mhSamples(J,:);
        chainResults.mhValue = chainResults.mhValue(J);
        if ifit==1
            [~,ir] = max(chainResults.mhValue);
        else
            ir = randi(size(chainResults.mhSamples,1));
        end
        x0 = exp(chainResults.mhSamples(ir,:));
        [parsMH,mlikelihoodMH] = Model.maximizeLikelihood(x0,fitOptions)
        if mlikelihoodMH<=mlikelihood
            pars = parsMH;
            disp('improved fit found')
        end

    catch
    end

    % Update Model and Make Plots of Results
    Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(pars);
end
Model.fspOptions.fspTol = 0.001;  % Set FSP error tolerance.
end

function [ModelZero,chainResults] = createDefaults(muLog10Prior)
for iPDO = 5:-1:1
    ModelZero{iPDO} = SSIT;    % Create SSIT instance and call it 'Model'.
    ModelZero{iPDO}.species = {'x1';'x2';'x3';'x4'}; % Set species names for bursting gene expression model:
    ModelZero{iPDO}.stoichiometry = [-1, 1, 0, 0, 0, 0;...
        1,-1,-1, 1, 0, 0;...
        0, 0, 1,-1, 0, 0;...
        0, 0, 0, 0, 1,-1];
    ModelZero{iPDO}.inputExpressions = {'ITrypt','(t<5)'};
    ModelZero{iPDO}.propensityFunctions = {'k12*x1';'k21*x2';'omega*x2';'1000*x3';'b*1000*x3*ITrypt';'g*x4'};
    ModelZero{iPDO}.initialCondition = [2;0;0;0];
    ModelZero{iPDO}.parameters = ({'k12',0.01;'k21',0.01;'omega',.2;'b',7.12;'g',0.012});
    ModelZero{iPDO}.parameters(:,2) = num2cell(10.^muLog10Prior);
    chainResults{iPDO}=[];
end
end
