% Benchmark Examples

Models = {'Toggle'};
clear benchmarks
for iM = 1:length(Models)
    Model = Generate_Model_from_Benchmark_Library(Models{iM});
    Model.propensityFilePrefix = Models{iM};
    benchmarks.(Models{iM}) = run_benchmarks(Model);
end

function Model = Generate_Model_from_Benchmark_Library(Name)
arguments 
    Name 
end

switch Name
    case 'Toggle'
        Model = SSIT;
        Model.parameters = {'kb',10;'ka',80;'M',20;'g',1};
        Model.species = {'LacI';'LamCI'};
        Model.stoichiometry = [1,-1,0, 0;...
            0, 0,1,-1];
        Model.propensityFunctions = {'kb+ka*M^3/(M^3+LamCI^3)';...
            'g*LacI';...
            'kb+ka*M^3/(M^3+LacI^3)';...
            'g*LamCI'};
        Model.initialCondition = [0;0];
        Model.customConstraintFuns = {'(LacI-3).^2.*(LamCI-3).^2'};
end

end

function benchmarks = run_benchmarks(Model,opts)
arguments
    Model
    opts.nSims = 100000;
end
% FSP solutions
Model.solutionScheme = 'fsp';
tic
Model = Model.formPropensitiesGeneral(Model.propensityFilePrefix);
benchmarks.writeFSPcodes = toc;

tic
[~,~,Model] = Model.solve;
benchmarks.initialFSPSolve = toc;

tic
[~,~,Model] = Model.solve;
benchmarks.subsequentFSPSolve = toc;

% SSA Solutions
Model.solutionScheme = 'ssa';
Model.ssaOptions.Nsims = 1;
tic
[~,~,Model] = Model.solve;
benchmarks.initialSSASolve_1run = toc;

Model.ssaOptions.Nsims = opts.nSims;
tic
[~,~,Model] = Model.solve;
benchmarks.(['subsequentSSASolve_',num2str(opts.nSims),'runs']) = toc;
opts.nSims;

%% ODE Solver
Model.solutionScheme = 'ode';
tic
[~,~,Model] = Model.solve;
benchmarks.i



end
