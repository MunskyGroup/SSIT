close all force
clear all
A = SSIT_GUI;  % Start up the app.

%% change to poisson model
A.ModelDropDown.Value = 'poisson_process.mat';
A.ModelDropDown.ValueChangedFcn(A,[]);

%% run SSA
A.PrintTimesEditField.Value = '[0 1 2 3 4 5]';
A.SsaRunButton.ButtonPushedFcn(A, [])
% you can ignore the error that pops up here for now.
%% save results
exportSsaData(A,'tmp_file.xlsx');

%% open SSA results in data fitting tab
fn = 'tmp_file.xlsx';
pn = [];
loadDataFromFile(A,fn,pn);
A.ParEstX1CheckBox.Value = 1;
A.ParEstX1DropDown.Value = {'x1'};
A.ParEstFitTimesList.Value = A.ParEstFitTimesList.Items;
userSelectCondition(A)

%% make plot of model and data
ssit.parest.updateModelSolveAndCompareToData(A);
makePlotOfData(A);

%% change to ham. mc.
A.FspErrorTolField.Value = inf;
% A.FittingAlgorithmDropDown.Value = 'Hamiltonian MC';
% A.FittingAlgorithmDropDown.Value = 'fminsearch';
% A.FittingAlgorithmDropDown.Value = 'Metropolis Hastings';
ssit.parest.fitModel2Data(A);
ssit.parest.updateModelSolveAndCompareToData(A);
makePlotOfData(A);