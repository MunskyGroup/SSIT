clear all
close all
nRounds = 8;
initialExperiment = [];
nFIMsamples = 10;

for rngSeed = 1

    %% Poisson Model With Chainging Input
    % truePars = {'kr',10;'gr',0.3;'kD',5};
    % inputLibrary = {{'IDex','1'},{'IDex','2'},{'IDex','3'},{'IDex','4'},{'IDex','5'},{'IDex','6'},...
    %     {'IDex','7'},{'IDex','8'},{'IDex','9'},{'IDex','10'}};
    % TestCases = struct('model',{'poisson','poisson'},...
    %                    'expDesign',{'fimopt','random'},...
    %                    'numCellsPerExperiment',{100,100},...
    %                    'figureSet',{1,1},...
    %                    'incrementAdd',{10,10},...
    %                    'truePars',{truePars,truePars},...
    %                    'saveFileNamePrefix',{'poisson_fim','poisson_rand'}, ...
    %                    'initialParameterGuess',{ones(1,3),ones(1,3)},...
    %                    'inputLibrary',{inputLibrary,inputLibrary});
    % 

    %% GR Model (Real)
    % truePars = {'kcn0',0.005;'kcn1',0.08;'knc',0.014;...
    %     'kg1',0.012;'gnuc',0.005;'MDex',10.44};
    load EricModel_MMDex
    truePars = {'kcn0',GRpars(1);'kcn1',GRpars(2);'knc',GRpars(4);...
        'kg1',GRpars(5);'gnuc',GRpars(7);'MDex',GRpars(8)};

    inputLibrary = {{'IDex','1'},{'IDex','10'},{'IDex','100'}};
    initialGuess = 10.^[-2 -1 -1 -2 -2 1];
    % initialGuess = [truePars{:,2}];
    TestCases = struct('model',{'dexgr','dexgr'},...
        'expDesign',{'fimopt','random'},...
        'numCellsPerExperiment',{80,80},...
        'figureSet',{1,1},...
        'incrementAdd',{10,10},...
        'truePars',{truePars,truePars},...
        'saveFileNamePrefix',{'gr_fim_real','gr_rand_real'}, ...
        'initialParameterGuess',{initialGuess,initialGuess},...
        'inputLibrary',{inputLibrary,inputLibrary},...
        'dataType',{'real','real'});
    initialExperiment = zeros(3,7);
    initialExperiment(3,[1,4,7]) = 20;
    initialExperiment(1:3,4) = 20;
    
    %% Burst Model
    % TestCases = struct('model',{'burst','burst'},...
    %     'expDesign',{'fimopt','random'},...
    %     'numCellsPerExperiment',{100,100}, ...
    %     'figureSet',{2,2},...
    %     'incrementAdd',{10,10},...
    %     'truePars',...
    %     {{'kon',0.1;'koff',0.2;'kr',10;'gr',0.3},{'kon',0.1;'koff',0.2;'kr',10;'gr',0.3}});

    %% Uncertain Mechanism Model
    % truePars1 = {'kon',0.1;'koff',0.2;'kr',10;'gr',0.3;'alph',1e-4};
    % truePars2 = {'kon',0.1;'koff',0.2;'kr',10;'gr',0.3;'alph',1e4};
    % alph = 0: Koff Control
    % alph  = inf: Kon Control
    % TestCases = struct('model',{'uncertainburst','uncertainburst','uncertainburst','uncertainburst'},...
    %     'expDesign',{'fimopt','random','fimopt','random'},...
    %     'numCellsPerExperiment',{100,100,100,100}, ...
    %     'figureSet',{1,1,2,2},...
    %     'incrementAdd',{20,20,20,20},...
    %     'truePars',{truePars1,truePars1,truePars2,truePars2},...
    %     'saveFileNamePrefix',{'koffBurstControlFIM','koffBurstControlRAND',...
    %                           'konBurstControlFIM','konBurstControlRAND'}, ...
    %     'initialParameterGuess',{ones(1,5),ones(1,5),ones(1,5),ones(1,5)});

    %% Both Model
    % TestCases = struct('model',{'poisson','poisson','burst','burst'},...
    %     'expDesign',{'fimopt','random','fimopt','random'},...
    %     'numCellsPerExperiment',{40,40,100,100}, ...
    %     'figureSet',{1,1,2,2},...
    %     'incrementAdd',{10,10,10,10},...
    %     'truePars',{{'kr',10;'gr',0.3},{'kr',10;'gr',0.3},...
    %     {'kon',0.1;'koff',0.2;'kr',10;'gr',0.3},{'kon',0.1;'koff',0.2;'kr',10;'gr',0.3}});

    %% Run Sequential Design Simulations
    for i = 1:length(TestCases)
        TestCases(i).saveName = iterativeExperimentRunnerMultiConditions(TestCases(i).model,...
            TestCases(i).dataType,...
            TestCases(i).expDesign,nRounds,rngSeed,rngSeed,...
            TestCases(i).incrementAdd,...
            TestCases(i).numCellsPerExperiment,...
            initialExperiment,nFIMsamples,...
            TestCases(i).truePars,...
            TestCases(i).saveFileNamePrefix,...
            TestCases(i).initialParameterGuess,...
            TestCases(i).inputLibrary);
    end
end

%% Load results and make plots
for i = 1:length(TestCases)
    load(TestCases(i).saveName)

    for round = 1:nRounds
        TestCases(i).detCov(round) = det(covLogMH{round});
        TestCases(i).pars(round,:) = parametersFound{round};
        
        TestCases(i).FIMpredNextExpt = FIMpredNextExpt;
        TestCases(i).MHResultsSaved = MHResultsSaved;

        TestCases(i).detFIMInv_Pred(round) = predictCov(FIMpredNextExpt{round});
        TestCases(i).detFIMTrueInv(round) = predictCov(FIMcurrentExptTrueSaved{round});

        TestCases(i).measurements = exptDesigns;
    end

    try
        clf(TestCases(i).figureSet*100+1)
        clf(TestCases(i).figureSet*100+2)
    catch
    end
end

cols = ['rkbg'];

%% Make figures of det(COV) vs experiment round
for j = 1:3
    for i = 1:length(TestCases)
        fset = TestCases(i).figureSet;
        f = figure(fset*100+1);
        set(f,'Name',[TestCases(i).model,' performance'])
        switch j
            case 1
                plot(1:nRounds,TestCases(i).detCov,'linewidth',2,'Color',cols(i)); hold on;
            case 2
                plot(2:nRounds,TestCases(i).detFIMInv_Pred(1:end-1),'--','linewidth',2,'Color',cols(i)); hold on;
            case 3
                plot(1:nRounds,TestCases(i).detFIMTrueInv,'-.','linewidth',2,'Color',cols(i)); hold on;
        end
        set(gca,"FontSize",16,'yscale','log')
        legend(TestCases([TestCases.figureSet]==fset).expDesign)

    end
end
xlabel('Experiment round');

%% Make Figures for MH Results and FIM Analyses
for i = 1:length(TestCases)
    fset = TestCases(i).figureSet;
    for iExpt = 2:length(FIMpredNextExpt)
        f = figure(fset*100+10+i*length(FIMpredNextExpt)+iExpt); clf;
        set(f,'Name',[TestCases(i).model,', ',TestCases(i).expDesign,', round ',num2str(iExpt)])
        SSIT.plotMHResultsStatic([],...
            TestCases(i).MHResultsSaved{iExpt}, ...
            TestCases(i).FIMpredNextExpt{iExpt-1},...
            'log','log10',f,false); 
        hold on;
        npars = size(TestCases(i).truePars,1);
        for iy = 1:npars-1
            for ix = iy+1:npars
                subplot(npars-1,npars-1,ix-1+(iy-1)*(npars-1))
                ylims = get(gca,'YLim');
                xlims = get(gca,'XLim');
                plot(log10(TestCases(i).truePars{ix,2})*[1,1],ylims,'k--')
                plot(xlims,log10(TestCases(i).truePars{iy,2})*[1,1],'k--')

            end
        end
    end
end

%% Make figure of parameter convergence
linspecs = ['--';'-.';'--';'-.'];
cols = ['rbrb'];
for i = 1:length(TestCases)
    fset = TestCases(i).figureSet;
    f = figure(fset*100+3);

    set(f,'Name',['Parameters, ' TestCases(i).model])
    semilogy(TestCases(i).pars,linspecs(i,:),'LineWidth',3,'Color',cols(i));hold on;
    for jpar = 1:size(TestCases(i).truePars,1)
        plot([1,size(TestCases(1).pars,1)],TestCases(i).truePars{jpar,2}*[1,1],'k--')
    end
    
    set(gca,"FontSize",16)
    xlabel('Experiment Iteration')
    ylabel('Parameter Estimates')

    set(gca,'ylim',[5e-2,2e1])   
end

%% Make figures of final experiment designs
for i = 1:length(TestCases)
    fset = TestCases(i).figureSet;
    f = figure(fset*100+2);
    set(f,'Name',[TestCases(i).model,' experiment design, ',TestCases(i).expDesign])
    measurements = 0;
    for j = 1:length(TestCases(i).measurements)
        measurements = measurements + TestCases(i).measurements{j};
    end
    bar(measurements); hold on;    
    set(gca,"FontSize",16)
    xlabel('Measurement Times')
    ylabel('Number of Cells Measured')
end

function predictedDetCov = predictCov(fimSet,parSet)
arguments
    fimSet
    parSet = [1:length(fimSet{1})];
end
detCovPreds = zeros(1,length(fimSet));
for i = 1:length(fimSet)
    detCovPreds(i) = det(fimSet{i}(parSet,parSet)^(-1));
end
predictedDetCov = mean(detCovPreds);
end

% plot(2:nRounds,detFIMTrueInv_FIM,'k',2:nRounds,detFIMTrueInv_FIM,'--k',2:nRounds,detFIMTrueInv_FIM,'-.k','linewidth',2)
% legend('|COV|','|FIM^-1|','|FIM^-1| True')
