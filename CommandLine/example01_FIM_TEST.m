%% BirthDeathTest
% In this script, we show how to set up and solve the FSP-FIM matrix with
% partial observations and probabilistic distortion.
clear all
ModelChoice = 'BirthDeath';  % Two species problem (mRNa and protein)
F = SSIT(ModelChoice);
F.solutionScheme = 'FSP';    % Set solution scheme to FSP.
F.fspOptions.fspTol = 1e-4;
[FSPsoln,bounds] = F.solve;  % Solve the FSP analysis
F.fspOptions.bounds = bounds;% Save bound for faster analyses 

%% FSP Solution Test
theoreticalMean = F.initialCondition*exp(-F.parameters{2,2}*F.tSpan) +...
    F.parameters{1,2}/F.parameters{2,2}*(1-exp(-F.parameters{2,2}*F.tSpan));

app.FspTabOutputs.solutions = FSPsoln.fsp;
app.FspPrintTimesField.Value = mat2str(F.tSpan);
solution = exportFSPResults(app);
err = abs(solution.Means - theoreticalMean')./theoreticalMean';
if max(err)<1e-4
    disp(['Mean(t) is correct within 0.01% at all t (',num2str(max(err),3),')'])
end

for t = 1:length(F.tSpan)
    fspPdf = solution.Marginals{t}{1}';
    theoreticalPdf = pdf('poiss',[0:length(fspPdf)-1],theoreticalMean(t));
    err(t) = sum(abs(theoreticalPdf-fspPdf));
end
if max(err)<0.001
     disp(['PDF is correct within one norm of 0.001 (',num2str(max(err),3),')'])
end

%% FSP Sensitivity Test
F.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = F.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

theoreticalMean = F.initialCondition*exp(-F.parameters{2,2}*F.tSpan) +...
    F.parameters{1,2}/F.parameters{2,2}*(1-exp(-F.parameters{2,2}*F.tSpan));
dk = F.parameters{1,2}*0.000001;
theoreticalMean2 = F.initialCondition*exp(-F.parameters{2,2}*F.tSpan) +...
    (F.parameters{1,2}+dk)/F.parameters{2,2}*(1-exp(-F.parameters{2,2}*F.tSpan));

for t = 1:length(F.tSpan)
    sensPdf = sensSoln.plotable.sensmdist{1,1,t}';
    theoreticalPdf = pdf('poiss',[0:length(sensPdf)-1],theoreticalMean(t));
    theoreticalPdf2 = pdf('poiss',[0:length(sensPdf)-1],theoreticalMean2(t));
    theoreticalSens = (theoreticalPdf2-theoreticalPdf)/dk;

    if max(theoreticalSens)>0
        err(t) = sum(abs(sensPdf-theoreticalSens)/max(abs(theoreticalSens)));
    else 
        err(t)=0;
    end
end
if max(err)<0.01
     disp(['PDF Sensitivity is correct within 1% of maximum (',num2str(max(err),3),')'])
end

%% Verify FIM using exact soln for Pure Birth Poisson Process.
F2 = SSIT('Empty');
F2.parameters = {'k',10};
F2.propensityFunctions = {'k'};
F2.stoichiometry = [1];
F2.species = {'x1'};
F2.initialCondition = 0;

F2.solutionScheme = 'FSP';    % Set solution scheme to FSP.
F2.fspOptions.fspTol = 1e-4;
[FSPsoln,F2.fspOptions.bounds] = F2.solve;  % Solve the FSP analysis

F2.solutionScheme = 'fspSens';    % Set solution scheme to FSP.
[sensSoln,bounds] = F2.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

F2.pdoOptions.unobservedSpecies = [];
F2.pdoOptions.PDO = [];
[fimResults] = F2.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.

% FIM = (dlam/dk)^2 * 1/lam = (d(k*t)/dk)^2 / (k*t) = t/k
analyticalFIM =  F2.tSpan/F2.parameters{1,2};

err = abs([fimResults{:}]-analyticalFIM)./analyticalFIM;
err(([fimResults{:}]-analyticalFIM)==0) = 0;

if max(err)<0.01
     disp(['FIM(t) is correct within 1% at all times (',num2str(max(err),3),')'])
end


%% Verify FIM for Binomial PDO (no dynamics)
F3 = SSIT('Empty');
F3.parameters = {};
F3.propensityFunctions = {'0'};
F3.stoichiometry = [1];
F3.species = {'x1'};
F3.initialCondition = 20;
F3.solutionScheme = 'FSP';    % Set solution scheme to FSP.
F3.fspOptions.fspTol = 1e-4;
[FSPsoln,F3.fspOptions.bounds] = F3.solve;  % Solve the FSP analysis

F3.solutionScheme = 'fspSens';    % Set solution scheme to FSP.
[sensSoln,bounds] = F3.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

pdoOptions.type = 'Binomial';
% Need to define loss parameter for each species S1, S2,...
pdoOptions.props.CaptureProbabilityS1 = 0.9; 
F3.pdoOptions.PDO = F3.generatePDO(pdoOptions,[],FSPsoln.fsp,true);
[fimResults_BinomialPDO] = F3.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.

% FIM = n / (p*(1-p)) % analytical form for binomial FI for 'p'.
analyticalFIM = F3.initialCondition/(pdoOptions.props.CaptureProbabilityS1*(1-pdoOptions.props.CaptureProbabilityS1));

err = abs([fimResults_BinomialPDO{:}]-analyticalFIM)/analyticalFIM;

if max(err)<0.01
     disp(['FIM(t) is correct within 1% at all times (',num2str(max(err),3),')'])
end

%% Verify FIM for Poisson PDO (no dynamics)
F3 = SSIT('Empty');
F3.parameters = {};
F3.propensityFunctions = {'0'};
F3.stoichiometry = [1];
F3.species = {'x1'};
F3.initialCondition = 20;
F3.solutionScheme = 'FSP';    % Set solution scheme to FSP.
F3.fspOptions.fspTol = 1e-4;
[FSPsoln,F3.fspOptions.bounds] = F3.solve;  % Solve the FSP analysis

F3.solutionScheme = 'fspSens';    % Set solution scheme to FSP.
[sensSoln,bounds] = F3.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

pdoOptions.type = 'Poisson';
% Need to define loss parameter for each species S1, S2,...
pdoOptions.props.NoiseMeanS1 = 10; 
F3.pdoOptions.PDO = F3.generatePDO(pdoOptions,[],FSPsoln.fsp,true);
[fimResults_PoissonPDO] = F3.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.

% FIM = 1 / lambda % analytical form for poisson FI for 'lambda'.
analyticalFIM = 1/pdoOptions.props.NoiseMeanS1;
 
err = abs([fimResults_PoissonPDO{:}]-analyticalFIM)/analyticalFIM;

if max(err)<0.01
      disp(['FIM(t) is correct within 1% at all times (',num2str(max(err),3),')'])
end

%% General verification of FIM
% Set up a model
F3 = SSIT('BirthDeath');
F3.tSpan = [0:10];
F3.solutionScheme = 'FSP';    % Set solution scheme to FSP.
F3.fspOptions.fspTol = 1e-4;
[FSPsoln,F3.fspOptions.bounds] = F3.solve;  % Solve the FSP analysis

% Generate some independent SSA trajectories:
F3.solutionScheme = 'SSA';    % Set solution scheme to SSA.
Nexp = 100;
Nt = length(F3.tSpan);
nSims = 100;
F3.ssaOptions.nSims = Nexp*nSims*Nt;
ssaSolns =  F3.solve;
%%
A = table;
for j=1:Nt
    A.t((j-1)*nSims+1:j*nSims) = F3.tSpan(j);
    for i = 1:Nexp
        for k=1:nSims
            for s = 1:size(ssaSolns.trajs,1)
                A.(['exp',num2str(i),'_s',num2str(s)])((j-1)*nSims+k) = ...
                    ssaSolns.trajs(s,j,(i-1)*Nt*nSims+(j-1)*nSims+k);
            end
        end
    end
end
writetable(A,'TEMP.csv')
    


