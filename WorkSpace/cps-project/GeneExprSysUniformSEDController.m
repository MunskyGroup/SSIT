% Gene Expression System Uniform SED Controller Class Definition
% This controller engages in uniform sequential experimental design, i.e.
% it distributes the remaining cells evenly among the remaining time
% points.
classdef GeneExprSysUniformSEDController < GeneExprSysSEDController
   methods(Access = ?GeneExprSysSEDController)
       function next_experiment = apply_ed_strategy(controller)
            % Uniform distribution of experiments
            matrixRows = controller.MaxRounds;
            matrixCols = controller.nInputs*controller.nT;
            uniformCell = floor(...
                ones(matrixRows, matrixCols) * ...
                controller.NumCellsPerExperiment/(matrixCols));
            for i=1:matrixRows
                if sum(uniformCell(i,:)) < controller.NumCellsPerExperiment
                    uniformCell(i,1:controller.NumCellsPerExperiment-sum(uniformCell(i,:)))=...
                        uniformCell(i,1:controller.NumCellsPerExperiment-sum(uniformCell(i,:)))+1;
                end
            end
            next_experiment = ...
                GeneExprSysCellDistribution(matrixRows, matrixCols);
            next_experiment.Matrix = uniform_cell;       
       end
    end 
end