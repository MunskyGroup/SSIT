%   'Determinant' - maximize the expected determinant of the FIM
%   'DetCovariance' - minimize the expected determinant of MLE covariance. 
%   'Smallest Eigenvalue' - maximize the smallest e.val of the
%       FIM
%   'Trace' - maximize the trace of the FIM
%   '[<i1>,<i2>,...]' - minimize the determinant of the inverse
%       FIM for the specified indices. All other parameters are
%       assumed to be free.
%   'TR[<i1>,<i2>,...]' - maximize the determinant of the  FIM
%       for the specified indices. Only the parameters in
%       obj.fittingOptions.modelVarsToFit are assumed to be
%       free.
classdef (Abstract) FIMOptimalityCriterion < matlab.mixin.Heterogeneous
    % Abstract classes define a template for a class, but you can't create
    % objects of this class. However you can define (non-abstract)
    % subclasses that inherit from it, and make objects of those classes.
    %
    % Inheriting from matlab.mixin.Heterogeneous means you can make arrays
    % of objects from the subclasses, with different subclasses in the same
    % array.

    methods(Abstract)
        metric_func = getMetricFunction(criterion)
    end
end