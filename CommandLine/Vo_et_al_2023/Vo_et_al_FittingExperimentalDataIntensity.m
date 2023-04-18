%% Vo_et_al_FittingExperimentalData.m
addpath('../')
TMP = SSIT(); % Load definition of the SSIT class for later use.
clear TMP

fileStr = 'FitResultsFeb15'; % name for saving the fitting results
vars.timeSet = [0,300]; %Data to fit: {'0','0_18','0_300','0_18_300'}
vars.doFit = 1;
vars.display = 'final';
vars.fitCases = [1:7];
vars.modelVarsToFit = [2:5]; % Variables allowed to change in fit.

%% Step 1 - Find the PDO parameters.
vars.doFit = 1;
vars.iter = 50000;
vars.pdoTimes = [0,300]; % Times to use in fitting the PDO
FittingFunctions(1,fileStr,vars);
%%
vars.doFit = 0;
results = FittingFunctions(1,fileStr,vars);
KSFISHspots = [(results.KS([1,5],1)+results.KS([1,5],3))/2,results.KS([1,5],2)]
KSMCPspots = [(results.KS([2,6],1)+results.KS([2,6],3))/2,results.KS([2,6],2)]
KSFISHintens = [(results.KS([7,3],1)+results.KS([7,3],3))/2,results.KS([7,3],2)]
KSMCPintens = [(results.KS([8,4],1)+results.KS([8,4],3))/2,results.KS([8,4],2)]

% KLDFISHspots = [(results.KLD([1,5],1)+results.KLD([1,5],3))/2,results.KLD([1,5],2)]
% KLDMCPspots = [(results.KLD([2,6],1)+results.KLD([2,6],3))/2,results.KLD([2,6],2)]
% KLDFISHintens = [(results.KLD([7,3],1)+results.KLD([7,3],3))/2,results.KLD([7,3],2)]
% KLDMCPintens = [(results.KLD([8,4],1)+results.KLD([8,4],3))/2,results.KLD([8,4],2)]

%% Step 2 - Fit the Data (fminsearch)
vars.doFit = 1;
vars.iter = 20; % Number of iterations in fminsearch
vars.display = 'final';  % Display type for fminsearch
vars.timeSet = [0,300]; % Data to fit: {'0','0_18','0_300','0_18_300'}
FittingFunctions(32,[fileStr,'_0_300'],vars); % Fitting the data at chosen time points
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
FittingFunctions(3,[fileStr,'_0_300'],vars)
%% Step 4 - Quantify Uncertainty (Met. Hast.)
vars.fitCases = [1,4,5,6,7];
vars.nMH = 1000;
vars.mhScaling = 0.3;
FittingFunctions(4,[fileStr,'_0_300'],vars)
% It is often helpful to run a few short chains (e.g., ~100) and adjust the proposal
% scaling until you get an acceptance rate of about 0.2-0.4.  This will
% lead to better mixing of the chain.
%% Step 5 - Compute FIM for subsequent experiment designs
vars.nFIMsamples = 1;
FittingFunctions(5,[fileStr,'_0_300'],vars)
