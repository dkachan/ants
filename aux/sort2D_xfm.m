function [connect2D] = sort2D_xfm(D,L,element_length)

global node

h = min(element_length)/2;

ny  = round(D/h+1);
nx  = round(L/h+1);

connect2D = zeros(ny,nx);

for j = 1:ny
    row  = find(abs(node(:,2)+D/2-h*(j-1))<1e-5);
    pos = round((node(row,1)+L/2+h)/h);
    
    connect2D(j,pos') = row';
end




