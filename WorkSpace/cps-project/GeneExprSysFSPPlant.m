% Gene Expression System FSP Plant Class Definition
classdef GeneExprSysFSPPlant < GeneExprSysPlant
    methods
        function data = RunNextExperiment(plant, experiment)
            arguments
                plant GeneExprSysPlant
                experiment GeneExprSysExperiment
            end
            curModel = plant.TrueModel; % Clone the model
            curModel.inputExpressions = experiment.Input;
            curModel = curModel.formPropensitiesGeneral(...
                [saveFileName,'_s',experiment.Input], true);
            fitOptions = optimset('Display','none','MaxIter',maxFitIter);
            curModel.fspOptions.fspTol = 1e-4;

            %% Generate Model Propensity Functions and Solve True Model
            [ModelSolution, curModel.fspOptions.bounds] = curModel.solve;
          
            %% FIM options
            fimScale = 'log'; % Maximize fim for log parameters
        
            %% True Model FIM
            curModel.fspOptions.fspTol = 1e-8;
            fimTrue = curModel.computeFIM([] ,fimScale);
        
            %% Verify that the true model and simulated data look correct.
            if showPlots
                dataFile = ['simData/TrueExperiment_',saveFileName,'_',num2str(inputIdx),'.csv'];
                nextExperiment = 100*ones(1,nT);
                curModel.ssaOptions.nSimsPerExpt = max(nextExperiment);
                curModel.ssaOptions.Nexp = 1;
                curModel.sampleDataFromFSP(curModel,dataFile)
                curModel = curModel.loadData(dataFile,dataToFit);
                curModel.makeFitPlot([],1);
            end

        end
    end
end