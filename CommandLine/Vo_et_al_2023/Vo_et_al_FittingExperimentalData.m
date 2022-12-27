addpath('../')
TMP = SSIT(); % Load definition of the SSIT class for later use.
clear TMP

fileStr = 'Name_for_fitting_case'; % name for savin ght efitting results
vars.timeSet = '0_300'; %Data to fit: {'0','0_18','0_300','0_18_300'}
vars.doFit = 1;
%% Step 1 - Find the PDOs.
vars.pdoTimes = [0,300]; % Times to use in fitting the PDO
FittingFunctionsCoLocalized(1,fileStr,vars)

%% Step 2 - Fit the Data (fminsearch)
vars.modelVarsToFit = [2:5]; % Variables allowed to change in fit.
vars.iter = 400; % Number of iterations in fminsearch
vars.display = 'final'; % Display type for fminsearch
switch vars.timeSet
    case '0'
        FittingFunctionsCoLocalized(2,fileStr,vars); %Fitting the data at t=0;
    case '0_18'
        FittingFunctionsCoLocalized(12,fileStr,vars); %Fitting the data at t=0,18;
    case '0_300'
        FittingFunctionsCoLocalized(22,fileStr,vars); %Fitting the data at t=0,300;
    case '0_18_300'
        FittingFunctionsCoLocalized(32,fileStr,vars); %Fitting the data at t=0,18,300;
end
% It is recommended that you run this a couple times before continuing to
% the next step.
%% Step 3 - Compute Fisher Information Matrix (current experiment design)
FittingFunctionsCoLocalized(3,fileStr,vars)
%% Step 4 - Quantify Uncertainty (Met. Hast.)
vars.nMH = 100;
FittingFunctionsCoLocalized(4,fileStr,vars)

%% Step 5 - Compute FIM for subsequent experiment designs
vars.nFIMsamples = 2;
FittingFunctionsCoLocalized(5,fileStr,vars)
