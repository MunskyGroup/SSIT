addpath('../CommandLine')

Poiss = SSIT;
Poiss.species = {'rna'};
Poiss.initialCondition = 0;
Poiss.propensityFunctions = {'kr';'gr*rna'};
Poiss.stoichiometry = [1,-1];
Poiss.parameters = ({'kr',10;'gr',1});
Poiss.tSpan = linspace(0,2,21);
Poiss.fspOptions.fspTol = 1e-5;
Poiss = Poiss.formPropensitiesGeneral('Poiss',true);

Poiss.solutionScheme = 'FSP';
[PoissSolution,Poiss.fspOptions.bounds] = Poiss.solve;

%%
clc
Poiss.sensOptions.solutionMethod = 'forward';

Poiss.solutionScheme = 'fspSens';
[PoissSolution,Poiss.fspOptions.bounds] = Poiss.solve;

Poiss.sensOptions.solutionMethod = 'finiteDifference';
SensSoln2 = Poiss.solve;

tic
Poiss.sensOptions.solutionMethod = 'forward';
Poiss.solutionScheme = 'fspSens';
[PoissSolution,Poiss.fspOptions.bounds] = Poiss.solve;
forwardTime = toc

tic
Poiss.sensOptions.solutionMethod = 'finiteDifference';
SensSoln2 = Poiss.solve;
finiteDiffTime = toc


[double(SensSoln2.sens.data{end}.S(1).data)-...
    double(PoissSolution.sens.data{end}.S(1).data)]