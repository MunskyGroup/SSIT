function app = changeFitType(app)
app.DataLoadingAndFittingTabOutputs.fitOptions = ssit.parest.propsStorage;

switch app.FittingAlgorithmDropDown.Value
    case 'fminsearch'
        app.DataLoadingAndFittingTabOutputs.fitOptions.props = ...
            optimset('display','iter',...
            'MaxIter',100);

    case 'simulated annealing'
        app.DataLoadingAndFittingTabOutputs.fitOptions.props = ssit.parest.anneal();
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.Verbosity = 2;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.Generator = @(x)(x+0.01*randn(size(x)));
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.InitTemp = 1;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.StopTemp = 1e-8;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.CoolSched = @(T)0.8*T;

    case 'genetic algorithm'
        app.DataLoadingAndFittingTabOutputs.fitOptions.props = optimoptions('ga');
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.UseParallel = true;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.Display = 'iter';
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.CrossoverFraction = 0.1;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.PopulationSize = 40;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.EliteCount = 1;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.MaxGenerations = 10;

    case 'particle swarm'
        app.DataLoadingAndFittingTabOutputs.fitOptions.props = optimoptions('particleSwarm');
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.UseParallel = true;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.Display = 'iter';
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.MaxIterations = 20;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.SwarmSize = 100;        
        
    case 'Metropolis Hastings'
        app.DataLoadingAndFittingTabOutputs.fitOptions.props = [];
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.isPropDistSymmetric=true;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.thin=1;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.numberOfSamples=1000;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.burnIn=100;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.progress=true;
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.proposalDistribution='@(x)x+0.01*randn(size(x))';
        app.DataLoadingAndFittingTabOutputs.fitOptions.props.numChains = 1;
        
end