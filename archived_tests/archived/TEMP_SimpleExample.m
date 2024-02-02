  %% Test Case 1 - a simple Poisson model
            testCase1.Poiss = SSIT;
            testCase1.Poiss.species = {'rna'};
            testCase1.Poiss.initialCondition = 0;
            testCase1.Poiss.propensityFunctions = {'kr';'gr*rna'};
            testCase1.Poiss.stoichiometry = [1,-1];
            testCase1.Poiss.parameters = ({'kr',10;'gr',1});
            testCase1.Poiss.tSpan = linspace(0,2,21);
            testCase1.Poiss.fspOptions.fspTol = 1e-5;
            testCase1.Poiss = testCase1.Poiss.formPropensitiesGeneral;

            [testCase1.PoissSolution,testCase1.Poiss.fspOptions.bounds] = testCase1.Poiss.solve;
            tic
            [testCase1.PoissSolution,testCase1.Poiss.fspOptions.bounds] = testCase1.Poiss.solve(testCase1.PoissSolution.stateSpace);
            testCase1.PoissSolution.time = toc;
