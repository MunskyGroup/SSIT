% Invariant Class Definition
classdef (Abstract) Invariant < matlab.mixin.Heterogeneous
    % Inheriting from matlab.mixin.Heterogeneous means you can make arrays
    % of objects from the subclasses, with different subclasses in the same
    % array.

    methods(Abstract)
        valid = predicate()
    end
end