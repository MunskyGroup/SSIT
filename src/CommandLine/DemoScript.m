%% Configurables and configurations

eic = ExperimentInputConfigurable();
eic.InputName = "Dex";
eic.Values = [0 2 4 6 8 10];
etc = ExperimentTimeConfigurable();
etc.Values = [0 30 60 90];
configs = multiplyConfigurations([eic etc]);

%% Designer specification

sed = SequentialExperimentDesigner();
sed.DataType = ExperimentalDataType.Empirical;
[~, dataFilename, ~] = fileparts(tempname(pwd));
sed.Model = DiscoverableModelFactory.createModel("GR", dataFilename);
sed.ExperimentalConfigurations = configs;

%% Round specification and design

round = SequentialExperimentRound(configs);
t = round.exportAsTable();
randomStrategy = RandomSEDStrategy();
round = randomStrategy.apportionObservations(round);
t = round.exportAsTable();
csv = round.exportAsCSV(tempname(pwd));

%% Import experimental data

opts = detectImportOptions(...
    [pwd '\ExampleDataSets\DUSP1_E_SSITcellresults_MG3_Abs6_Jun11_mg_abs.csv']);
opts.VariableTypes = strrep(opts.VariableTypes, 'char', 'string');
data = readtable([pwd ...
    '\ExampleDataSets\DUSP1_E_SSITcellresults_MG3_Abs6_Jun11_mg_abs.csv'], ...
    opts);