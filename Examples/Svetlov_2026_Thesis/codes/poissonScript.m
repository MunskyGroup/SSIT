%% Ground-truth model

ModelTrue = MBExperimentModel;
ModelTrue.species = {'rna'};
ModelTrue.initialCondition = 0;
ModelTrue.propensityFunctions = {'kr*IDex/(kD+IDex)'; 'gr*rna'};
ModelTrue.stoichiometry = [1, -1];
ModelTrue.parameters = {'kr', 10; 'gr', 0.3; 'kD', 5};

muLog10Prior = [0,0,0];
sigLog10Prior = [1 1 1];
ModelTrue = ModelTrue.setPriorsLog10(muLog10Prior, sigLog10Prior);
%dataToFit = {'rna','exp1_s1'};
%fitParameters = [1:3];


%% Guessed model

ModelGuess = ModelTrue;
ModelGuess.parameters = {'kr', 1; 'gr', 1; 'kD', 1};

%% Configurations

inputConfigurable = MBExperimentInputConfigurable();
inputConfigurable.InputName = "IDex";
inputConfigurable.Values = 1:10;

timeConfigurable = MBExperimentTimeConfigurable();
nT = 21;
timeConfigurable.Values = linspace(0, 10, nT);

configs = MBExperimentDesigner.multiplyConfigurables(...
    [inputConfigurable timeConfigurable]);

% if isempty(inputLibrary)
%     ModelTrue.inputExpressions = {'IDex','5'};
%     ModelTrue = ModelTrue.formPropensitiesGeneral([datFileName,'_S',num2str(ind)],true);
%     ModelTrue = {ModelTrue};
% else
%     TMP = cell(1,length(inputLibrary));
%     for iInput = 1:length(inputLibrary)
%         TMP{iInput} = ModelTrue;
%         TMP{iInput}.inputExpressions = inputLibrary{iInput};
%         TMP{iInput} = TMP{iInput}.formPropensitiesGeneral([datFileName,'_s',num2str(iInput)],true);
%     end
%     ModelTrue = TMP;
% end

%% Design Strategy

strategy = MBExperimentFIMOptimizedDesignStrategy();
info = strategy.FIMInfo();
info.NumberOfObservationsPerExperiment = 60;
%fimMetric = 'D-cov'; %'Determinant'';

%% Initial Design

design = MBExperimentDesign(configs);
[Obs, nonTimeConfigurations, times] = design.getAsObservationMatrix();
Obs(1, [1, 7, 14, 21]) = round(info.NumberOfObservationsPerExperiment / 4);
design = design.setFromObservationMatrix(Obs, nonTimeConfigurations, times);

%% Designer Creation and Operation

designer = MBExperimentDesigner(...
    [inputConfigurable timeConfigurable], ...
    ModelGuess, ...
    "", ... % idString
    strategy, ...
    design, ...
    ModelTrue, ...
    1);
designer.NumberOfMHSamplesForProduction = 1000;
designer.UseEmpiricalData = false;
designer.NumberOfFIMSamples = 10;
%nSamplesMH = 1000; % Number of MH Samples to run

%nRounds = 8;
%nFIMsamples = 10;
%datType = 'simulated';

designer.performNextRound();
design = designer.designNextRound();