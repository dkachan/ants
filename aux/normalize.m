function out = normalize( arr )
%NORMALIZE Summary of this function goes here
%   Detailed explanation goes here

out = arr./(sqrt(sum(arr.^2, 2))*ones(1,size(arr,2)));

%{
[A, B] = size(arr); 
C = zeros(A,1); 
for i=1:B
    C = C + arr(:,i).^2;
end

out = arr./(C.^(1/2)*ones(1,B));
%}