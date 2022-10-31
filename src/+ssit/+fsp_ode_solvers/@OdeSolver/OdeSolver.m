classdef (Abstract) OdeSolver 
    %ODESOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)        
    end
    
    methods
        function obj = OdeSolver()
        %ODESOLVER Construct an instance of this class                
        end
    end
    
    methods (Abstract)
        [tExport, solutionsNow, fspStopStatus] = solve(obj, tStart, tOut, initSolution, rhs, jac, fspErrorCondition)                       
    end
end

