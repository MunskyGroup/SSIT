% Gene Expression System Random SED Controller Class Definition
% This controller engages in random sequential experimental design, i.e.
% it distributes the remaining cells randomly among the remaining time
% points.
classdef GeneExprSysRandomSEDController < GeneExprSysSEDController
   methods(Access = ?GeneExprSysSEDController)
       function next_experiment = apply_ed_strategy(controller)
            % Random distribution of experiments
            % Numbers of cells must be integral
            N = controller.NumCellsPerExperiment/controller.cellBatchSize;
            matrixRows = controller.MaxRounds;
            matrixCols = controller.nInputs*controller.nT;
            randomCell = GeneExprSysCellDistribution(...
                matrixRows, matrixCols);
            for i = 1:matrixRows
                n = randi(matrixCols,1,N);
                for j = n
                    randomCell.Matrix(i,j) = ...
                        randomCell.Matrix(i,j) + controller.cellBatchSize;
                end
            end
            next_experiment.Matrix = randomCell;
       end
    end 
end