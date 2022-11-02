clear all
close all

ModelChoice = 'BirthDeath';    % One species problem
% ModelChoice = 'ToggleSwitch';  % Two species problem (non-linear toggle switch)
% ModelChoice = 'CentralDogmaTV';  % Two species problem (mRNa and protein) with time varying transcription rate
% ModelChoice = 'Repressilator'; % Three species problem 
% ModelChoice = 'BurstingSpatialCentralDogma';   % Four species problem (gene state, nuclear mRNa, cytoplasmic mRNA, and protein)

%% Create SSIT Model
F2 = SSIT(ModelChoice);

%% Extending a model:
% ModelChoice = 'BurstingSpatialCentralDogma';  
% F2 = SSIT(ModelChoice);
% F2 = F2.addSpecies('x5');
% F2 = F2.addParameter({'kpact',0.2;'gpact',0.1});
% F2 = F2.addReaction('kpact*x4',[0;0;0;-1;1]); 
% F2 = F2.addReaction('gpact*x5',[0;0;0;0;-1]); 
% % Note -- this 5-species reaction that will require an FSP
% % solution size of  >10^7 ODEs.  It will likely take several minutes
% % to solve on most computers.

%% Solve using FSP
F2.tSpan = [-1:1:10];
F2.initialTime = -1;
F2.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
F2.fspOptions.verbose=1;
[FSPsoln,bounds] = F2.solve;  % Solve the FSP analysis

%% Solve again using the bounds from the last solution
% If we start with the bounds computed in the first analysis, the solution is
% often much faster.
F2.fspOptions.bounds = bounds;% Save bound for faster analyses 
[FSPsoln] = F2.solve(FSPsoln.stateSpace);  % Solve the FSP analysis

%% Make plots of FSP solution
F2.makePlot(FSPsoln,'meansAndDevs') % Make some plots.

%% Solve Sensitivity using FSP
F2.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = F2.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

%% Make plot of sensitivities
F2.makePlot(sensSoln.plotable,'marginals') % Plot marginal sensitivities

%% Run SSA
F2.solutionScheme = 'SSA';
SSASoln = F2.solve;
