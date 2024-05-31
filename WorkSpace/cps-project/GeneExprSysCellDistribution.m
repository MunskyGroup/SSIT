% Gene Expression System Experiment Cell Distribution Class Definition
classdef GeneExprSysCellDistribution
    properties
        Matrix (:,:) integer {mustBeNonnegative}       
    end
    methods
        function cd = GeneExprSysCellDistribution(rows, cols)
            arguments
                rows (1,1) integer {mustBePositive}
                cols (1,1) integer {mustBePositive}
            end
            cd.Matrix = zeros(rows, cols);
        end
        
        function cd_out = plus(cd1, cd2)
            if ~(size(cd1.Matrix) == size(cd2.DistributionMatrix))
                error('Cannot add cell distribution matrices of unequal dimensions!')                
            end
            cd_out = GeneExprSysCellDistribution;
            cd_out.DistributionMatrix = ...
                cd1.Matrix + cd2.DistributionMatrix;
        end
    end
end