classdef miscelaneousTests < matlab.unittest.TestCase
    % this tests:
    %       1) Creation of models using SBML.
    properties
        Model
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(tc)
            addpath('../CommandLine')
            tc.Model = SSIT;
         end  
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)

        function loadModelFromSBML(tc)
            % Tests the loading of a model from SBML.
            tc.Model = tc.Model.createModelFromSBML('../SBML_test_cases/00010/00010-sbml-l1v2.xml');
            [fspSoln] = tc.Model.solve;
            tc.Model.makePlot(fspSoln,'meansAndDevs')
        end

        function testNonlinearLogicTimeVarying(tc)

                        %% Test Case 3 - 2 Species Poisson Model
            TwoDNonLinearTV = SSIT;
            TwoDNonLinearTV.species = {'rna1','rna2'};
            TwoDNonLinearTV.initialCondition = [0;0];
            TwoDNonLinearTV.propensityFunctions = {'kr1';'gr1*rna1*(1/(1+(rna2/M)^eta*Ir))';'kr2';'gr2*rna2'};
            TwoDNonLinearTV.stoichiometry = [1,-1,0,0;0,0,1,-1];
            TwoDNonLinearTV.parameters = ({'kr1',20;'gr1',1;'kr2',18;'gr2',1;'M',20;'eta',5});
            TwoDNonLinearTV.inputExpressions = {'Ir','t>1'};
            TwoDNonLinearTV.tSpan = linspace(0,2,21);
            TwoDNonLinearTV = TwoDNonLinearTV.formPropensitiesGeneral('TwoNonLinTV');

            [TwoDNonLinearTVSolution,TwoDNonLinearTV.fspOptions.bounds] = TwoDNonLinearTV.solve;
            tic
            [TwoDNonLinearTVSolution,TwoDNonLinearTV.fspOptions.bounds] = TwoDNonLinearTV.solve(TwoDNonLinearTVSolution.stateSpace);
            TwoDNonLinearTVSolution.time = toc;

            % %% ODE model for Two Poisson process
            % TwoDNonLinearTVODE = TwoDNonLinearTV;
            % TwoDNonLinearTVODE.solutionScheme = 'ode';
            % TwoDNonLinearTVODE = TwoDNonLinearTVODE.formPropensitiesGeneral('TwoDNonLinearTVODE');  


        end

    end
end