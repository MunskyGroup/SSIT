% Gene Expression System SED Controller Class Definition
% This controller engages in sequential experimental design according to
% a particular design strategy, e.g. random distribution of cells among
% remaining time points.
classdef (Abstract) GeneExprSysSEDController < GeneExprSysController
    properties (SetAccess = private)
        cellBatchSize (1,1) integer {mustBePositive} = 10
        cellDistributions (1,1) GeneExprSysCellDistribution
        covMH (1,:) cell
        covLogMH (1,:) cell
        covLogFIM_Prediction (1,:) cell
        covFIM_Prediction (1,:) cell
        exptDesigns (1,:) cell
        FIMpredNextExpt (1,:) cell
        FIMcurrentExptSaved (1,:) cell
        FIMcurrentExptTrueSaved (1,:) cell
        FIMsamples (1,1) integer {mustBePositive} = 10
        MHResultsSaved (1,:) cell
        nInputs (1,1) integer {mustBePositive}
        nT (1,1) integer {mustBePositive}
        parametersFound (1,:) cell
        showPlots (1,1) boolean = false
    end

    methods
        % Need constructor to allow for extra properties
        function c = GeneExprSysSEDController(...
                model, rounds, batchSize, showPlots)
            superargs = {model, rounds};
            % Call superclass constructor
            c@GeneExprSysController(superargs{:});

            if (nargin > 2)
                c.cellBatchSize = batchSize;
            end
            if (nargin > 3)
                c.showPlots = showPlots;
            end

            % Initialize cell arrays for saving results
            c.covMH = cell(1,rounds);
            c.covLogMH = cell(1,rounds);
            c.covLogFIM_Prediction = cell(1,rounds);
            c.covFIM_Prediction = cell(1,rounds);
            c.parametersFound = cell(1,rounds);
            c.exptDesigns = cell(1,rounds);
            c.FIMpredNextExpt = cell(1,rounds);
            c.FIMcurrentExptSaved = cell(1,rounds);
            c.FIMcurrentExptTrueSaved = cell(1,rounds);
            c.MHResultsSaved = cell(1,rounds);

            % Initialize cell distribution matrix for possible experiments
            c.nT = length(c.EstimatedModel.tSpan);
            c.nInputs = length(c.InputLibrary);
            c.cellDistributions = GeneExprSysCellDistribution(...
                c.nInputs, c.nT);
            c.cellDistributions.Matrix(...
                1, [1, round(nT/3), round(2*nT/3), nT]) = ...
                round(c.NumCellsPerExperiment/4);
        end
    end

    methods(Sealed, Access = ?GeneExprSysController)
        function experiment = SelectNextExperiment(controller)
            nextRound = controller.CompletedRounds + 1;
            clc
            disp(['Round: ',num2str(nextRound)])

            % Adjust the time span. All previous data are fixed (they
            % will have been provided by the Plant, and thus are
            % correct). We can only control the signal and cell
            % distribution for the remaining time interval, which we
            % will reset to a zero origin.
            newTimeSpan = ...
                controller.EstimatedModel.tSpan(nextRound:controller.MaxRounds);
            if controller.CompletedRounds > 0
                newTimeSpan = newTimeSpan - controller.EstimatedModel.tSpan(...
                    controller.CompletedRounds);
            end

            %% Create model variants for possible inputs
            
            % Generate simulated data
            
                

            % 
            
            for iInput = 1:controller.nInputs
                
            end

            
            

            
        end % SelectNextExperiment

        function updateModel(controller, data)
        end

        function objFun = getObjective(x,Mods,hasData,stateSpace)
        % This function takes a group of models and common parameter set 
        % and then calculates the log-likelihood function.
            arguments
                x
                Mods
                hasData
                stateSpace = [];           
            end
            if isempty(stateSpace)
                stateSpace = cell(1,length(Mods));
            end
            objFun = 0;
            includePrior = true;
            for i = 1:length(Mods)
                % Check if there is data associated with that model.
                if hasData(i)
                    % We only want to include the prior in the first likelihood
                    % calculation.  Otherwise it will be amplified again in every data
                    % set.
                    if includePrior
                        objFun = objFun + ...
                            Mods{i}.computeLikelihood(exp(x),stateSpace{i});
                        includePrior = false;
                    else
                        TMP = Mods{i};
                        TMP.fittingOptions.logPrior = [];
                        objFun = objFun + ...
                            TMP.computeLikelihood(exp(x),stateSpace{i});            
                    end
                end
            end
        end % getObjective
    end % Sealed functions

    methods(Abstract, Access = ?GeneExprSysSEDController)
        next_experiment = apply_ed_strategy(controller)
    end
end
            

