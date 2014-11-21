classdef Obstacle < matlab.mixin.Heterogeneous
    %OBSTACLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type;
    end
    
    methods
        function obj = Obstacle(type)
            obj.type = type;
        end
    end
    
end

