classdef propsStorage < handle  %  <===== THIS LINE MAKES IT A HANDLE CLASS
    properties 
        % We can store any amount of data in this empty structure
        props = struct(); 
    end
    
    methods 
        function obj = propsStorage(obj)
            return;
        end
    end
end