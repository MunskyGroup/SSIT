% Gene Expression System FIM Controller Class Definition
% This calculates the "uncertainty volume" (the determinant of the
% covariance matrix) for the posterior distribution after each experiment
% and selects the next experiment such that the expected subsequent
% uncertainty volume is minimized.
classdef GeneExprSysFIMController < GeneExprSysController
    properties (SetAccess = private)
        covMH (1,:) cell
        covLogMH (1,:) cell
        covLogFIM_Prediction (1,:) cell
        covFIM_Prediction (1,:) cell
        exptDesigns (1,:) cell
        FIMpredNextExpt (1,:) cell
        FIMcurrentExptSaved (1,:) cell
        FIMcurrentExptTrueSaved (1,:) cell
        FIMsamples (1,1) integer {mustBePositive} = 10
        MHResultsSaved (1,:) cell
        parametersFound (1,:) cell
    end

    methods
        % Need constructor to allow for extra properties
        function c = GeneExprSysFIMController(model, rounds, FIMsamples)
            superargs = {model, rounds};
            % Call superclass constructor
            c@GeneExprSysController(superargs{:});

            c.FIMsamples = FIMsamples;

            % Initialize cell arrays for saving results
            c.covMH = cell(1,rounds);
            c.covLogMH = cell(1,rounds);
            c.covLogFIM_Prediction = cell(1,rounds);
            c.covFIM_Prediction = cell(1,rounds);
            c.parametersFound = cell(1,rounds);
            c.exptDesigns = cell(1,rounds);
            c.FIMpredNextExpt = cell(1,rounds);
            c.FIMcurrentExptSaved = cell(1,rounds);
            c.FIMcurrentExptTrueSaved = cell(1,rounds);
            c.MHResultsSaved = cell(1,rounds);
        end
    end

    methods(Access = ?GeneExprSysController)
        function experiment = SelectNextExperiment(controller)
        end
        function updateModel(controller, data)
        end
    end
end
            
function [TestCases,finalExperimentDesign] = sequentialExptDesignBatchRunner(iModel,iDesign,rngSeed, makePlots, testing, plotTimes)

switch iDesign
    case 1
        design = 'fimopt';
    case 2
        design = 'random';
end
incrementAdd = 10;
switch iModel
    case 1 % Poisson
        %% Poisson Model With Chainging Input
        truePars = {'kr',10;'gr',0.3;'kD',5};
        inputLibrary = {{'IDex','1'},{'IDex','2'},{'IDex','3'},{'IDex','4'},{'IDex','5'},{'IDex','6'},...
            {'IDex','7'},{'IDex','8'},{'IDex','9'},{'IDex','10'}};
        model = 'poisson';
        numCellsPerExperiment = 60;
        initialParGuess = ones(1,3);
        datType = 'simulated';
        nInputs = length(inputLibrary);
        nT = 21;
        initialExperiment = zeros(nInputs,nT);
        initialExperiment(1,[1,7,14,21]) = round(numCellsPerExperiment/4);

    case 2 % GR Model (Simulated Data)
        %% GR Model (Simulated Data)
        load MMDexPars truePars

        inputLibrary = {{'IDex','1'},{'IDex','10'},{'IDex','100'}};
        initialParGuess = [truePars{:,2}];
        % initialParGuess = 10.^[-2 -1 -2 -2 -2 1];
        model = 'dexgr';
        numCellsPerExperiment = 200;
        datType = 'simulated';

        initialExperiment = zeros(3,6);
        initialExperiment(3,[1:6]) = 50;
        initialExperiment(1:3,6) = 50;
end


TestCases = struct('model',{model},...
    'expDesign',{design},...
    'numCellsPerExperiment',{numCellsPerExperiment},...
    'figureSet',{iModel},...
    'incrementAdd',{incrementAdd},...
    'truePars',{truePars},...
    'saveFileNamePrefix',{['mod_',num2str(iModel),model,'_',num2str(iDesign),design,'_',num2str(rngSeed),'seed']}, ...
    'initialParameterGuess',{initialParGuess},...
    'inputLibrary',{inputLibrary},...
    'dataType',{datType});


%% Run Sequential Design Simulations
TestCases.saveName = iterativeExperimentRunnerMultiConditions(TestCases.model,...
        TestCases.dataType,...
        TestCases.expDesign,nRounds,rngSeed,jobID,...
        TestCases.incrementAdd,...
        TestCases.numCellsPerExperiment,...
        initialExperiment,nFIMsamples,...
        TestCases.truePars,...
        TestCases.saveFileNamePrefix,...
        TestCases.initialParameterGuess,...
        TestCases.inputLibrary,...
        testing);

if ~makePlots
    return
end


if iModel==4||iModel==5
    parset = [1:5];
elseif iModel==8
    parset = [2,5,6];
else
    parset = [1:size(truePars,1)];
end


%% Load results and make plots
load(TestCases.saveName)
nRounds = 0;
while nRounds<length(covLogMH)&&~isempty(covLogMH{nRounds+1})
    nRounds=nRounds+1;
end

for iRound = 1:nRounds
    TestCases.detCov(iRound) = det(covLogMH{iRound}(parset,parset));
    TestCases.pars(iRound,:) = parametersFound{iRound};

    TestCases.FIMpredNextExpt = FIMpredNextExpt;
    TestCases.MHResultsSaved = MHResultsSaved;

    [TestCases.detFIMInv_Pred(iRound),...
        TestCases.stdDetFIMInv_Pred(iRound),...
        TestCases.minDetFIMInv_Pred(iRound),...
        TestCases.maxDetFIMInv_Pred(iRound)] = ...
        predictCov(FIMpredNextExpt{iRound},parset);
    
    [TestCases.detFIMInv_Post(iRound),...
        TestCases.stdDetFIMInv_Post(iRound),...
        TestCases.minDetFIMInv_Post(iRound),...
        TestCases.maxDetFIMInv_Post(iRound)] = ...
        predictCov(FIMcurrentExptSaved{iRound},parset);

    TestCases.detFIMTrueInv(iRound) = predictCov(FIMcurrentExptTrueSaved{iRound},parset);
    
    TestCases.FIMTrue(iRound) = FIMcurrentExptTrueSaved{iRound};

    TestCases.measurements = exptDesigns;
end

cols = ['rkbg'];
%% Make figures of det(COV) vs experiment round
for j = 1:3
    fset = TestCases.figureSet;
    f = figure(fset*100+1);
    set(f,'Name',[TestCases.model,' performance'])
    switch j
        case 1
            plot(1:nRounds,TestCases.detCov,'linewidth',2,'Color',cols(iDesign)); hold on;
        case 2
            plot(2:nRounds,TestCases.detFIMInv_Pred(1:end-1),'--','linewidth',2,'Color',cols(iDesign)); hold on;
            % errorbar(2:nRounds,TestCases.detFIMInv_Pred(1:end-1),...
            %     TestCases.stdDetFIMInv_Pred(1:end-1),...
            %     '--','linewidth',2,'Color',cols(iDesign)); hold on;
            % errorbar(2:nRounds,TestCases.detFIMInv_Pred(1:end-1),...
            %     TestCases.detFIMInv_Pred(1:end-1) - TestCases.minDetFIMInv_Pred(1:end-1),...
            %     TestCases.maxDetFIMInv_Pred(1:end-1) - TestCases.detFIMInv_Pred(1:end-1),...
            %     '--','linewidth',2,'Color',cols(iDesign)); hold on;
        case 3
            if strcmp(datType,'simulated')
                plot(1:nRounds,TestCases.detFIMTrueInv,'s','linewidth',2,'MarkerSize',16,'Color',cols(iDesign),'MarkerFaceColor',cols(iDesign)); hold on;
            else
                plot(1:nRounds,TestCases.detFIMInv_Post,'s','linewidth',2,'MarkerSize',16,'Color',cols(iDesign),'MarkerFaceColor',cols(iDesign)); hold on;
            end          
    end
    set(gca,"FontSize",16,'yscale','log')
    legend(TestCases([TestCases.figureSet]==fset).expDesign)

end
xlabel('Experiment round');

%% Make Figures for MH Results and FIM Analyses
fset = TestCases.figureSet;
for iRound = 1:nRounds
    f = figure(fset*100+10+iDesign*length(FIMpredNextExpt)+iRound); 
    set(f,'Name',[TestCases.model,', ',TestCases.expDesign,', round ',num2str(iRound)])
    if iRound==1
        SSIT.plotMHResultsStatic([],...
            TestCases.MHResultsSaved{iRound}, ...
            TestCases.FIMTrue{iRound},...
            'log','log10',f,true);
    else
        SSIT.plotMHResultsStatic([],...
            TestCases.MHResultsSaved{iRound}, ...
            TestCases.FIMpredNextExpt{iRound-1},...
            'log','log10',f,true);
    end
end

%% Make figure of parameter convergence
linspecs = ['--';'-.';'--';'-.'];
cols = ['rbrb'];
fset = TestCases.figureSet;
f = figure(fset*100+3);

set(f,'Name',['Parameters, ' TestCases.model])
semilogy(TestCases.pars,linspecs(iDesign,:),'LineWidth',3,'Color',cols(iDesign));hold on;
for jpar = 1:size(TestCases.truePars,1)
    plot([1,size(TestCases.pars,1)],TestCases.truePars{jpar,2}*[1,1],'k--')
end


%% Make figures of final experiment designs
% fset = TestCases.figureSet;
% f = figure(fset*100+2);
% set(f,'Name',[TestCases.model,' experiment design, ',TestCases.expDesign])
finalExperimentDesign = initialExperiment;
if iModel==8
    for j = 1:min(nRounds-1,5)
        finalExperimentDesign = finalExperimentDesign + TestCases.measurements{j};
    end

else

    for j = 1:nRounds-1
        finalExperimentDesign = finalExperimentDesign + TestCases.measurements{j};
    end
end
% bar(finalExperimentDesign'); hold on;
% set(gca,"FontSize",16)
% xlabel('Measurement Times')
% ylabel('Number of Cells Measured')
end

function [predictedDetCov,stdvDetCov,minDetCov,maxDetCov] = predictCov(fimSet,parSet)
arguments
    fimSet
    parSet = [1:length(fimSet{1})];
end
detCovPreds = zeros(1,length(fimSet));
for i = 1:length(fimSet)
    covEst = fimSet{i}^(-1);
    detCovPreds(i) = det(covEst(parSet,parSet));
end
predictedDetCov = mean(detCovPreds);
stdvDetCov = std(detCovPreds);
minDetCov = min(detCovPreds);
maxDetCov = max(detCovPreds);
end

% plot(2:nRounds,detFIMTrueInv_FIM,'k',2:nRounds,detFIMTrueInv_FIM,'--k',2:nRounds,detFIMTrueInv_FIM,'-.k','linewidth',2)
% legend('|COV|','|FIM^-1|','|FIM^-1| True')