addpath(genpath('../../src'));


for iModel = 1:3
    switch iModel
        case 1 % Poisson TV model (size 10)
            disp('Model 1')
            Model = SSIT;
            Model.species = {'rna'};
            Model.initialCondition = 0;
            Model.stoichiometry = [1,-1];
            Model.parameters = ({'kr',10;'gr',1});
            Model.tSpan = linspace(0,2,21);
            Model.fspOptions.fspTol = 1e-5;

            Model.propensityFunctions = {'kr*Ig';'gr*rna'};
            Model.inputExpressions = {'Ig','t>1'};
            Model = Model.formPropensitiesGeneral('PoissTV');

        case 2 % Poisson TV model (size 1000)
            disp('Model 2')
            Model = SSIT;
            Model.species = {'rna'};
            Model.initialCondition = 0;
            Model.stoichiometry = [1,-1];
            Model.parameters = ({'kr',1000;'gr',1});
            Model.tSpan = linspace(0,2,21);
            Model.fspOptions.fspTol = 1e-5;

            Model.propensityFunctions = {'kr*Ig';'gr*rna'};
            Model.inputExpressions = {'Ig','t>1'};
            Model = Model.formPropensitiesGeneral('PoissTV');

        case 3 % Poisson TV model (size 100000)
            disp('Model 3')
            Model = SSIT;
            Model.species = {'rna'};
            Model.initialCondition = 0;
            Model.stoichiometry = [1,-1];
            Model.parameters = ({'kr',100000;'gr',1});
            Model.tSpan = linspace(0,2,21);
            Model.fspOptions.fspTol = 1e-5;

            Model.propensityFunctions = {'kr*Ig';'gr*rna'};
            Model.inputExpressions = {'Ig','t>1'};
            Model = Model.formPropensitiesGeneral('PoissTV');

    end

    [ModelSolution,Model.fspOptions.bounds] = Model.solve;

    %%
% Initialize a variable to store the case number
selectedCase = 0;

    for jSolnScheme = 1:6
        switch jSolnScheme
            case 1
                selectedCase = 1;
                disp('ode45')
                Model.fspOptions.odeSolutionScheme = 'ode45';
            case 2
                selectedCase = 2;
                disp('ode23s')
                Model.fspOptions.odeSolutionScheme = 'ode23s';
            case 3
                selectedCase = 3;
                disp('ode23s with Sparse Pattern')
                Model.fspOptions.odeSolutionScheme = 'ode23s_SparseWithJPattern';
            case 4
                selectedCase = 4;
                disp('ode23t')
                Model.fspOptions.odeSolutionScheme = 'ode23t';
            case 5
                selectedCase = 5;
                disp('ode15s')
                Model.fspOptions.odeSolutionScheme = 'ode15s';
            case 6
                selectedCase = 6;
                disp('ode113')
                Model.fspOptions.odeSolutionScheme = 'ode113';
        end

        tic
        for i=1:20
            [ModelSolution,Model.fspOptions.bounds] = Model.solve(ModelSolution.stateSpace);
        end
        time = toc;
        disp(['The computation time for case ' num2str(selectedCase) ' is: ' num2str(time) ' seconds'])

    end
end