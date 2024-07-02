%% Model Settings for SSIT Command Line Run

%% Import data into SSIT Model.
fileName = 'EricDataJan23_2024/pdoCalibrationData_EricIntensity_DexSweeps.csv';
linkedSpecies = {'rna','RNA_DUSP1_nuc'};
conditions = {'Dex_Conc','100'};
obj = obj.loadData(fileName, linkedSpecies, conditions);

executeRoutine = @(SSIT)executionRoutine(SSIT);

%% Functions to be executed.
function outputs = executionRoutine(ssitModel)

%% Run the fit
fitOptions = optimset('Display','none','MaxIter',300);
fitOptions.suppressFSPExpansion = true;
  
indsPars = ssitModel.fittingOptions.modelVarsToFit;
pars = [ssitModel.parameters{indsPars,2}];
pars = ssitModel.maximizeLikelihood(pars, fitOptions);
ssitModel.parameters(indsPars,2) = num2cell(pars);

%% Run MH Search
%TODO

%% Compute FIM
%TODO

%% Find Optimal Next Experiment
%TODO

%% Save Results in Updated Model File
%%TODO


outputs.model = ssitModel;
outputs.finished = true;
end