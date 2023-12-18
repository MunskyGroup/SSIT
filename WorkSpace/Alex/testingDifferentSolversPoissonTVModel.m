addpath(genpath('../../src'));


for iModel = 1:3
    switch iModel
        case 1 % Poisson TV model
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

        case 2  

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

        case 3
            disp('Model 3')
            Model = SSIT;
            Model.species = {'rna'};
            Model.initialCondition = 0;
            Model.stoichiometry = [1,-1];
            Model.parameters = ({'kr',10000;'gr',1});
            Model.tSpan = linspace(0,2,21);
            Model.fspOptions.fspTol = 1e-5;

            Model.propensityFunctions = {'kr*Ig';'gr*rna'};
            Model.inputExpressions = {'Ig','t>1'};
            Model = Model.formPropensitiesGeneral('PoissTV');


    end

    [ModelSolution,Model.fspOptions.bounds] = Model.solve;

    %%
    for jSolnScheme = 1:3
        switch jSolnScheme
            case 1
                disp('ode45')
                Model.fspOptions.odeSolutionScheme = 'ode45';
            case 2
                disp('ode23s')
                Model.fspOptions.odeSolutionScheme = 'ode23s';
            case 3
                disp('ode23s with Sparse Pattern')
                Model.fspOptions.odeSolutionScheme = 'ode23s_SparseWithJPattern';
        end

        tic
        for i=1:20
            [ModelSolution,Model.fspOptions.bounds] = Model.solve(ModelSolution.stateSpace);
        end
        time = toc

    end
end

