%% Model Settings for SSIT Command Line Run

%% Import data into SSIT Model.
% fileName = 'EricDataJan23_2024/pdoCalibrationData_EricIntensity_DexSweeps.csv';
linkedSpecies = {'rna','RNA_DUSP1_nuc';
    'rna','RNA_DUSP1_nuc';
    'rna','RNA_DUSP1_nuc';
    'rna','RNA_DUSP1_nuc'};
conditions = {'Dex_Conc','100';
    'Dex_Conc','10';
    'Dex_Conc','1';
    'Dex_Conc','0.3'};
freeParameters = {[1:4];[1:4];[1:4];[1:4]};

Model = obj.loadData(fileName, linkedSpecies, conditions);
Pars = [Model.parameters{Model.fittingOptions.modelVarsToFit,2}];

%% Remove prior from the model -- since this should only be applied once.
logPrior = @(x)Model.fittingOptions.logPrior(log(x));
Model.fittingOptions.logPrior = [];

%%  Load and Associate with DUSP1 smFISH Data
ModelGroup = cell(size(linkedSpecies,1),1);
boundGuesses = cell(size(linkedSpecies,1),1);
ModelGroupParameterMap = freeParameters;

for i = 1:size(linkedSpecies,1)
    ModelGroup{i} = Model.loadData(fileName,...
        linkedSpecies(i,:),...
        conditions(i,:)); 

    % Set Dex concentration in model
    ModelGroup{i}.parameters{13,2} = str2num(conditions{i,2});
        
    boundGuesses{i} = Model.fspOptions.bounds;
    % ModelGroup{i} = ModelGroup{i}.formPropensitiesGeneral(['Model_',num2str(i),'_FSP']);
end

combinedModel = SSITMultiModel(ModelGroup,ModelGroupParameterMap,logPrior);
combinedModel = combinedModel.initializeStateSpaces(boundGuesses);
combinedModel = combinedModel.updateModels(Pars,false);

executeRoutine = @(SSIT)executionRoutine(SSIT);


%% Functions to be executed.
function outputs = executionRoutine(combinedModel)

fitOptions = optimset('Display','iter','MaxIter',300);
Pars = combinedModel.parameters;

%% Run the fit
Pars = combinedModel.maximizeLikelihood(...
    Pars, fitOptions);
combinedModel = combinedModel.updateModels(Pars,false);

%% Run MH Search
%TODO

%% Compute FIM
%TODO

%% Find Optimal Next Experiment
%TODO

%% Save Results in Updated Model File
%%TODO

outputs.model = combinedModel;
outputs.finished = true;
end