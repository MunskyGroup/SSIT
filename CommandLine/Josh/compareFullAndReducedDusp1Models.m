%% Full Model (EGRNT)
tmp = load('complex_dusp1_model');
EGRNT = tmp.Model;

EGRNT.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
EGRNT.fspOptions.fspTol = 1e-8;  % Set FSP error tolerance.
[EGRNTsoln,EGRNT.fspOptions.bounds] = EGRNT.solve;  % Solve the FSP analysis
[EGRNTsoln,EGRNT.fspOptions.bounds] = EGRNT.solve;  % Solve the FSP analysis

EGRNT.ssaOptions.nSimsPerExpt = 100;
EGRNT.ssaOptions.Nexp = 200; 
EGRNT.sampleDataFromFSP(EGRNTsoln,'full_dusp1_model_testC.csv'); 

EGRNT.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
EGRNT = EGRNT.loadData('full_dusp1_model_testC.csv',{'x3','exp1_s3'});
EGRNT.initialTime = 0;
EGRNT.fittingOptions.timesToFit = ones(1,length(EGRNT.tSpan),'logical');
EGRNT.makeFitPlot

%% Reduced model (SGRS)
clc
SGRS = SSIT;
SGRS.species = {'x1';'x2'};  % GRnuc, geneOn, dusp1
SGRS.initialCondition = [0;0];
SGRS.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
SGRS.inputExpressions = {'IGR','kcn0/knc+(t>=0)*kcn1/(r1-knc)*(exp(-knc*t)-exp(-r1*t))'};
SGRS.stoichiometry = [1,-1,0,0;0,0,1,-1];
SGRS.parameters = EGRNT.parameters;
SGRS.tSpan = EGRNT.tSpan;
SGRS.fspOptions.initApproxSS = true;

SGRS.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
SGRS.fspOptions.fspTol = 1e-8;  % Set FSP error tolerance.
[SGRSsoln,SGRS.fspOptions.bounds] = SGRS.solve;  % Solve the FSP analysis
[SGRSsoln,SGRS.fspOptions.bounds] = SGRS.solve;  % Solve the FSP analysis
%% Plot Comparison of Full and Reduce Model
SGRS.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
SGRS = SGRS.loadData('full_dusp1_model_testC.csv',{'x2','exp1_s3'});
SGRS.initialTime = 0;
SGRS.fittingOptions.timesToFit = ones(1,length(SGRS.tSpan),'logical');
SGRS.makeFitPlot