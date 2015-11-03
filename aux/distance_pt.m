function [ d  ] = distance_pt( x,y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
d = ((x(1)-y(:,1)).^2+(x(2)-y(:,2)).^2).^(1/2);

end

