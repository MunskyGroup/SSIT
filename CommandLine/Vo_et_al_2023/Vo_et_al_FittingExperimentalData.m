% Vo_et_al_FittingExperimentalData.m
% addpath('../')
TMP = SSIT(); % Load definition of the SSIT class for later use.
clear TMP

fileStr = 'Name_for_fitting_case'; % name for saving the fitting results
vars.timeSet = [0,300]; %Data to fit: {'0','0_18','0_300','0_18_300'}
vars.doFit = 1;
%% Step 1 - Find the PDO parameters.
vars.pdoTimes = [0,300]; % Times to use in fitting the PDO
FittingFunctionsCoLocalized(1,fileStr,vars)

%% Step 2 - Fit the Data (fminsearch)
vars.modelVarsToFit = [2:5]; % Variables allowed to change in fit.
vars.iter = 400; % Number of iterations in fminsearch
vars.display = 'iter';  % Display type for fminsearch
vars.timeSet = [0,300]; % Data to fit: {'0','0_18','0_300','0_18_300'}
FittingFunctionsCoLocalized(32,fileStr,vars); % Fitting the data at chosen time points
% It is pretty common that you will get stuck in a local minimum during
% this fit.  To help with that, we run a few fits using the different PDOs
% and exchanging initial parametr guesses. 
% It is recommended that you run this a couple times before continuing to
% the next step. Also, it is recommended to redo this a couple times after
% completing Step 4 in cases  a better parameter set was found using the
% MH. To check how well the fit is doing, you can run the first few cells
% of the code "Vo_et_al_PlottingResults" (make sure to change the file
% names to your current fits).
%% Step 3 - Compute Fisher Information Matrix (current experiment design)
FittingFunctionsCoLocalized(3,fileStr,vars)
%% Step 4 - Quantify Uncertainty (Met. Hast.)
vars.nMH = 3000;
vars.mhScaling = 0.6;
FittingFunctionsCoLocalized(4,fileStr,vars)
% It is often helpful to run a few short chains (e.g., ~100) and adjust the proposal
% scaling until you get an acceptance rate of about 0.2-0.4.  This will
% lead to better mixing of the chain.
%% Step 5 - Compute FIM for subsequent experiment designs
vars.nFIMsamples = 4;
FittingFunctionsCoLocalized(5,fileStr,vars)
