close all
clear 
clc
%% Specifying parameters
nx=100;                           %Number of steps in space(x)
ny=100;                           %Number of steps in space(y)       
nt=5000;                           %Number of time steps 
dt=0.001;                         %Width of each time step
dx=2/(nx-1);                     %Width of space step(x)
dy=2/(ny-1);                     %Width of space step(y)
x=0:dx:2;                        %Range of x(0,2) and specifying the grid points
y=0:dy:2;                        %Range of y(0,2) and specifying the grid points
ants_0=zeros(nx,ny);                  %Preallocating u
ants_1=zeros(nx,ny);                  %Preallocating u
Fx=zeros(nx,ny);                  %Preallocating Fx
Fy=zeros(nx,ny);                  %Preallocating Fy
FObsx=zeros(nx,ny);                  %Preallocating Fx
FObsy=zeros(nx,ny);                  %Preallocating Fy
Back2Nx=zeros(nx,ny);                  %Preallocating Fx
Back2Ny=zeros(nx,ny);                  %Preallocating Fy
Obs=zeros(nx,ny);
un=zeros(nx,ny);                 %Preallocating un
zeta = 1;
vis=0.1;                         %Diffusion coefficient/viscocity
UW=0;                            %x=0 Dirichlet B.C 
UE=0;                            %x=L Dirichlet B.C 
US=0;                            %y=0 Dirichlet B.C 
UN=0;                            %y=L Dirichlet B.C 
UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/dn=UnN)
food_x = 1.5;
food_y = 1.5;
nest_x = 0.25;
nest_y = 0.25;




%% Initial Conditions
for i=1:nx
    for j=1:ny
        if ((nest_y-0.1<=y(j))&&(y(j)<=nest_y+0.1)&&(nest_x-0.1<=x(i))&&(x(i)<=nest_x+0.1))
            ants_0(i,j)=2;
        else
            ants_0(i,j)=0;
        end
        Fx(i,j) = 0;
        Fy(i,j) = 0;
        n = 2*[(x(i)-nest_x), (y(j)-nest_y)]/norm([(x(i)-nest_x), (y(j)-nest_y)]);
        Back2Nx(i,j) = -4*n(1);
        Back2Ny(i,j) = -4*n(2);
    end
end
Back2Nx(isnan(Back2Nx)) = 0;
Back2Ny(isnan(Back2Ny)) = 0;

%% Pheromones Initial Conditions
figure;
hold on
xlim([0 2])
ylim([0 2])
i = 1;
[X, Y] = meshgrid(x,y);
while (1)
    [xin(i), yin(i)] = ginput (1) ;
    if xin(i)>0;
        plot (xin(i), yin(i), 'o');
    end
    if xin(i)<0;
        xin(i)=[];
        yin(i)=[];
        break; end;
    i = i + 1;
end
curvexy(:,1)=xin;
curvexy(:,2)=yin;

mapxy = [X(:) Y(:)];
[xy,distance_to_trail_map_vec,~] = distance2curve(curvexy,mapxy,'Spline');

U = -1*reshape(distance_to_trail_map_vec,[nx,ny]).^(1/1.5);

% % Obstacles
% figure;
% hold on
% xlim([0 2])
% ylim([0 2])
% i = 1;
% [X, Y] = meshgrid(x,y);
% while (1)
%     [xin(i), yin(i)] = ginput (1) ;
%     if xin(i)>0;
%         plot (xin(i), yin(i), 'o');
%     end
%     if xin(i)<0;
%         xin(i)=[];
%         yin(i)=[];
%         break; end;
%     i = i + 1;
% end
% curvexy(:,1)=xin;
% curvexy(:,2)=yin;
% 
% mapxy = [X(:) Y(:)];
% [xy,distance_to_trail_map_vec,~] = distance2curve(curvexy,mapxy,'linear');
% 
% Obs = -1*reshape(distance_to_trail_map_vec,[nx,ny]);


%% Calculating the field variable for each time step



i=2:nx-1;
j=2:ny-1;

FObsx(i,j) = (U(i+1,j)-U(i-1,j))/(2*dx);
FObsy(i,j) = (U(i,j+1)-U(i,j-1))/(2*dy);
% 
% 
% figure
% h=surf(x,y,U','EdgeColor','none');       %plotting the field variable
% shading interp
% axis ([0 2 0 2 -2 2])

figure(1)
hold on
plot_arc(0,2*pi,nest_y,nest_x,0.1)
plot_arc(0,2*pi,food_y,food_x,0.1)

for it=0:nt
    ants_0n=ants_0;
    ants_1n=ants_1;
    
%     subplot(1,2,1)
    h=surf(x,y,ants_0,'EdgeColor','none');
    h2=surf(x,y,10*ants_1,'EdgeColor','none');
%     [C,h3] = contour(x,y,ants_0',[0.35,0.35]);
    axis ([0 2 0 2])
    caxis([0 0.45])
    drawnow;
    
    
    
%     subplot(1,2,2)
%     h2=surf(x,y,ants_1','EdgeColor','none');
%     axis ([0 2 0 2])
%     caxis([0 1])
%     drawnow;
    
    
    source = ants_0;
    source(abs(x-food_x)>0.05,:)=0;
    source(:,abs(y-food_y)>0.05)=0;
    
    source2 = ants_1;
    source2(abs(x-nest_x)>0.05,:)=0;
    source2(:,abs(y-nest_y)>0.05)=0;
    loading_rate = 0;
    unloading_rate = 0;
    pheromones = ants_1.^(1/4);
    
    Fx(i,j) = 2*(pheromones(i+1,j)-pheromones(i-1,j))/(2*dx);
    Fy(i,j) = 2*(pheromones(i,j+1)-pheromones(i,j-1))/(2*dy);
    
    ants_0(i,j)=ants_0n(i,j)+(vis*dt*(ants_0n(i+1,j)-2*ants_0n(i,j)+ants_0n(i-1,j))/(dx*dx))+(vis*dt*(ants_0n(i,j+1)-2*ants_0n(i,j)+ants_0n(i,j-1))/(dy*dy))...
        - dt*zeta^(-1)*(FObsx(i,j).*(ants_0n(i+1,j)-ants_0n(i-1,j))/(2*dx)  +  ants_0n(i,j).*(FObsx(i+1,j)-FObsx(i-1,j))/(2*dx))...
        - dt*zeta^(-1)*(FObsy(i,j).*(ants_0n(i,j+1)-ants_0n(i,j-1))/(2*dy)  +  ants_0n(i,j).*(FObsy(i,j+1)-FObsy(i,j-1))/(2*dy))...
        - dt*zeta^(-1)*(Fx(i,j).*(ants_0n(i+1,j)-ants_0n(i-1,j))/(2*dx)  +  ants_0n(i,j).*(Fx(i+1,j)-Fx(i-1,j))/(2*dx))...
        - dt*zeta^(-1)*(Fy(i,j).*(ants_0n(i,j+1)-ants_0n(i,j-1))/(2*dy)  +  ants_0n(i,j).*(Fy(i,j+1)-Fy(i,j-1))/(2*dy)) - loading_rate*source(i,j) + unloading_rate*source2(i,j);
    

    ants_1(i,j)=ants_1n(i,j)+(vis*dt*(ants_1n(i+1,j)-2*ants_1n(i,j)+ants_1n(i-1,j))/(dx*dx))+(vis*dt*(ants_1n(i,j+1)-2*ants_1n(i,j)+ants_1n(i,j-1))/(dy*dy)) ...
        - dt*zeta^(-1)*(FObsx(i,j).*(ants_1n(i+1,j)-ants_1n(i-1,j))/(2*dx)  +  ants_1n(i,j).*(FObsx(i+1,j)-FObsx(i-1,j))/(2*dx))...
        - dt*zeta^(-1)*(FObsy(i,j).*(ants_1n(i,j+1)-ants_1n(i,j-1))/(2*dy)  +  ants_1n(i,j).*(FObsy(i,j+1)-FObsy(i,j-1))/(2*dy))...
        - dt*zeta^(-1)*(Back2Nx(i,j).*(ants_1n(i+1,j)-ants_1n(i-1,j))/(2*dx)  +  ants_1n(i,j).*(Back2Nx(i+1,j)-Back2Nx(i-1,j))/(2*dx))...
        - dt*zeta^(-1)*(Back2Ny(i,j).*(ants_1n(i,j+1)-ants_1n(i,j-1))/(2*dy)  +  ants_1n(i,j).*(Back2Ny(i,j+1)-Back2Ny(i,j-1))/(2*dy))+ loading_rate*source(i,j) - unloading_rate*source2(i,j);
    
    
%     Boundary conditions
%     Dirichlet:
%     ants_0(1,:)=UW;
%     ants_0(nx,:)=UE;
%     ants_0(:,1)=US;
%     ants_0(:,ny)=UN;
%     
%     ants_0(x-0.8>0,:) = 0;
%     
%     ants_1(1,:)=UW;
%     ants_1(nx,:)=UE;
%     ants_1(:,1)=US;
%     ants_1(:,ny)=UN;
%     Neumann:
    ants_0(1,:)=ants_0(2,:)-UnW*dx;
    ants_0(nx,:)=ants_0(nx-1,:)+UnE*dx;
    ants_0(:,1)=ants_0(:,2)-UnS*dy;
    ants_0(:,ny)=ants_0(:,ny-1)+UnN*dy;
    
    ants_1(1,:)=ants_1(2,:)-UnW*dx;
    ants_1(nx,:)=ants_1(nx-1,:)+UnE*dx;
    ants_1(:,1)=ants_1(:,2)-UnS*dy;
    ants_1(:,ny)=ants_1(:,ny-1)+UnN*dy;
%     }
if it<nt
    delete(h);
    delete(h2);
end
% delete(h3);

end