function plot_mesh_dk( X, connect ,style)
%PLOT_MESH_DK Summary of this function goes here
%   Detailed explanation goes here

hold on
ord =[1 2 3 4 1];
xpt = zeros(size(connect,1),length(ord));
ypt = zeros(size(connect,1),length(ord));
for e=1:size(connect,1)
    for n=1:length(ord)
      xpt(e, n)=X(connect(e,ord(n)),1);
      ypt(e, n)=X(connect(e,ord(n)),2);      
    end
end
plot(xpt',ypt',style)
axis equal

end

