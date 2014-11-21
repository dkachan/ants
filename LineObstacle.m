classdef LineObstacle < Obstacle
    %LINEOBSTACLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pt1;
        pt2;
        width;
    end
    
    methods
        function obj = LineObstacle(pt1,pt2,width,type)
            obj@Obstacle(type);
            obj.width = width;
            obj.pt1 = pt1;
            obj.pt2 = pt2;
        end
    end
    
end

