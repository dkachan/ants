classdef CircleObstacle < Obstacle
    %CIRCLEOBSTACLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        center;
        radius;
    end
    
    methods
        function obj = CircleObstacle(center,radius)
            obj@Obstacle(ObstacleTypes.Circle);
            obj.center = center;
            obj.radius = radius;
        end
    end
    
end

