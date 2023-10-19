clear all
close all

%% (1) Creating and Modifying SSIT Model

%%      (1A) Choosing a Pre-Made Model
ModelChoice = 'BirthDeath';    % One species problem
% ModelChoice = 'ToggleSwitch';  % Two species problem (non-linear toggle switch)
% ModelChoice = 'CentralDogmaTV';  % Two species problem (mRNa and protein) with time varying transcription rate
% ModelChoice = 'Repressilator'; % Three species problem 
% ModelChoice = 'BurstingSpatialCentralDogma';   % Four species problem (gene state, nuclear mRNa, cytoplasmic mRNA, and protein)
F2 = SSIT(ModelChoice);

%%      (1B) Extending a model with additional reactions
% ModelChoice = 'BurstingSpatialCentralDogma';  
% F2 = SSIT(ModelChoice);
% F2 = F2.addSpecies('x5');
% F2 = F2.addParameter({'kpact',0.2;'gpact',0.1});
% F2 = F2.addReaction('kpact*x4',[0;0;0;-1;1]); 
% F2 = F2.addReaction('gpact*x5',[0;0;0;0;-1]); 
% % Note -- this 5-species reaction that will require an FSP
% % solution size of  >10^7 ODEs.  It will likely take several minutes
% % to solve on most computers.

%% (2) Solving a Model
%%      (2A) Solve using FSP
F2 = F2.formPropensitiesGeneral('BasicModel');
F2.tSpan = [-1:1:10];
F2.initialTime = -1;
F2.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
F2.fspOptions.verbose=1;
[FSPsoln,F2.fspOptions.bounds] = F2.solve;  % Solve the FSP analysis

%%          (2A.1) Solve FSP Model again using the bounds from the last solution
% If we start with the bounds computed in the first analysis, the solution is
% often much faster.
[FSPsoln] = F2.solve(FSPsoln.stateSpace);  % Solve the FSP analysis

%%          (2A.2) Make plots of FSP solution
F2.makePlot(FSPsoln,'meansAndDevs',[],[],1,{'linewidth',3,'color',[0,1,1]}) % Make plot of mean vs. time.
F2.makePlot(FSPsoln,'marginals',[],[],2,{'linewidth',3,'color',[0,0,1]}) % Make plot of mean vs. time.

%%      (2B) Solve using SSA
F3 = F2;
F3.solutionScheme = 'SSA';
SSASoln = F3.solve;

%%          (2B.1) Make plots of SSA solution
F3.makePlot(SSASoln,'trajectories',[],[],4) % Make some plots.
F2.makePlot(FSPsoln,'meansAndDevs',[],[],4,...
    {'linewidth',4,'color',[0,1,1],'Marker','s','MarkerSize',20}) % Add FSP Solution to plot.

%% (3) Sensitivity Analysis using FSP
F4 = F3;
F4.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = F4.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

%%      (3A) Make plot of sensitivities
F4.makePlot(sensSoln,'marginals',[],[],4,{'linewidth',3,'color',[0,0,1]}) % Plot marginal sensitivities

