close all
clear
clc

addpath('../aux');
%% Specifying parameters
nx=100;                           %Number of steps in space(x)
ny=100;                           %Number of steps in space(y)
nt=100000;                           %Number of time steps
dt=0.00001;                         %Width of each time step
dx=2/(nx-1);                     %Width of space step(x)
dy=2/(ny-1);                     %Width of space step(y)
x=0:dx:2;                        %Range of x(0,2) and specifying the grid points
y=0:dy:2;                        %Range of y(0,2) and specifying the grid points
ants_0=zeros(nx,ny);                  %Preallocating u
ants_1=zeros(nx,ny);                  %Preallocating u
Fx_0=zeros(nx,ny);                  %Preallocating Fx_0
Fy_0=zeros(nx,ny);                  %Preallocating Fy_0
FObsx=zeros(nx,ny);                  %Preallocating Fx
FObsy=zeros(nx,ny);                  %Preallocating Fy_0
Fx_1=zeros(nx,ny);                  %Preallocating Fx
Fy_1=zeros(nx,ny);                  %Preallocating Fy_0
distance_to_trail_map_vec_tot = zeros(nx*ny,1);  
Obs=zeros(nx,ny);
U=ones(nx,ny);
un=zeros(nx,ny);                 %Preallocating un
zeta_0 = 1;                      %mobility
zeta_1 =1;                       %mobility
vis_0=1;                       %Diffusion coefficient/viscocity
vis_1=1;                       %Diffusion coefficient/viscocity
UW=0;                            %x=0 Dirichlet B.C
UE=0;                            %x=L Dirichlet B.C
US=0;                            %y=0 Dirichlet B.C
UN=0;                            %y=L Dirichlet B.C
UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/dn=UnN)
loading_rate = 1;         %self explanatory
unloading_rate = 1;       %self explanatory
interac_exponent = 1  ;           %self explanatory
food_x = 1.2; %y
food_y = 1; %x
nest_x = 0.8; %y
nest_y = 1; %x
maxF0 = 1e20;                    %Force bound
maxF1 = 1e20;                     %Force bound
alpha0 = 1;
alpha1 = 1;


%% Initial Conditions
for i=1:nx
    for j=1:ny
        if ((nest_y-0.1<=y(j))&&(y(j)<=nest_y+0.1)&&(nest_x-0.1<=x(i))&&(x(i)<=nest_x+0.1))
            ants_0(i,j)=2;
        else
            ants_0(i,j)=0;
        end
        Fx_0(i,j) = 0;
        Fy_0(i,j) = 0;
        %         n = 2*[(x(i)-nest_x), (y(j)-nest_y)]/norm([(x(i)-nest_x), (y(j)-nest_y)]);
        Fx_1(i,j) = 0;
        Fy_1(i,j) = 0;
    end
end
[X, Y] = meshgrid(x,y);
mu_nest = [nest_y nest_x];
mu_food = [food_y food_x];
Sigma = [.003 0; 0 0.003];
G_nest = mvnpdf([X(:) Y(:)],mu_nest,Sigma);
G_food = mvnpdf([X(:) Y(:)],mu_food,Sigma);
nest_centered_gaussian = reshape(G_nest,length(x),length(y));
food_centered_gaussian = reshape(G_food,length(x),length(y));
nest_centered_gaussian = nest_centered_gaussian/max(max(nest_centered_gaussian));
food_centered_gaussian = food_centered_gaussian/max(max(food_centered_gaussian));

nest_centered_gaussian(nest_centered_gaussian<1e-2) = 0;
food_centered_gaussian(food_centered_gaussian<1e-2) = 0;

ants_0 = 1*reshape(nest_centered_gaussian,length(x),length(y));

% Fx_1(isnan(Fx_1)) = 0;
% Fy_1(isnan(Fy_1)) = 0;

%% Pheromones Initial Conditions
% figure;
% hold on
% xlim([0 2])
% ylim([0 2])

[X, Y] = meshgrid(x,y);
mapxy = [X(:) Y(:)];
obstaclesx =[];
obstaclesy = [];
% while 1
%     i = 1;
%     while (1)
%         [xin(i), yin(i)] = ginput (1) ;
%         if xin(i)>0;
%             plot (xin(i), yin(i), 'o');
%         end
%         if xin(i)<0;
%             xin(i)=[];
%             yin(i)=[];
%             break; end;
%         i = i + 1;
%     end
%     curvexy(:,1)=xin;
%     curvexy(:,2)=yin;
%     
%     [xy,distance_to_trail_map_vec,~] = distance2curve(curvexy,mapxy,'linear');
%     
%     distance_to_trail_map_vec_tot = [distance_to_trail_map_vec_tot  distance_to_trail_map_vec];
%     obstaclesx = [obstaclesx xin'];
%     obstaclesy = [obstaclesy yin'];
%     [xi, yi] = ginput (1) ;
%     if xi > 2;
%         break
%     end
% end
% 
% for i = 1:nx*ny
%     distance_to_trail_map_vec(i,1) = min(distance_to_trail_map_vec_tot(i,2:end));
% end
% 
% U = reshape(distance_to_trail_map_vec,[nx,ny]).^(1/1.5);
% save('U2.mat','U','obstaclesx','obstaclesy')
% load('U2.mat')

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
% Obs = 1*reshape(distance_to_trail_map_vec,[nx,ny]);

create_movie = true;

if create_movie
    mov2 = VideoWriter('ants_smoluchowski3.avi');
    open(mov2);
end

%% Calculating the field variable for each time step



i=2:nx-1;
j=2:ny-1;

FObsx(i,j) = -(Obs(i+1,j)-Obs(i-1,j))/(2*dx);
FObsy(i,j) = -(Obs(i,j+1)-Obs(i,j-1))/(2*dy);
%
%
% figure
% h=surf(x,y,U','EdgeColor','none');       %plotting the field variable
% shading interp
% axis ([0 2 0 2 -2 2])

hFig = figure('Color',[1 1 1]);
set(hFig, 'Position', [100 50 1000 400])
hold on
plot_arc(0,2*pi,nest_y,nest_x,0.1)
plot_arc(0,2*pi,food_y,food_x,0.1)


for it=0:nt
    ants_0n=ants_0;
    ants_1n=ants_1;
    
    if mod(it,100) == 0
        subplot(1,2,1)
        hold on
        h=surf(x,y,ants_0,'EdgeColor','none');
        plot3(obstaclesx,obstaclesy,2*ones(size(obstaclesx)),'r-')
        axis ([0 2 0 2])
        caxis([0 2])
        drawnow;
        subplot(1,2,2)
        hold on
        h2=surf(x,y,ants_1,'EdgeColor','none');
        plot3(obstaclesx,obstaclesy,2*ones(size(obstaclesx)),'r-')
        axis ([0 2 0 2])
%         caxis([0 0.0004])
        drawnow;
        if create_movie
            currFrame = getframe(hFig);
            writeVideo(mov2,currFrame);
        end
        pause(0.01)
    end
    caxis([0 0.4])
    
    ants_1(ants_1<0) = 0;
    ants_0(ants_0<0) = 0;
    
    pheromones_1 = ants_1.^(interac_exponent);
    Fx_0(i,j) = alpha0*(pheromones_1(i+1,j)-pheromones_1(i-1,j))/(2*dx);
    Fy_0(i,j) = alpha0*(pheromones_1(i,j+1)-pheromones_1(i,j-1))/(2*dy);
    
    Fx_0(Fx_0>maxF0)=maxF0;
    Fy_0(Fy_0>maxF0)=maxF0;
    Fx_0(Fx_0<-maxF0)=-maxF0;
    Fy_0(Fy_0<-maxF0)=-maxF0;
    
    pheromones_0 = ants_0.^(interac_exponent);
    Fx_1(i,j) = alpha1*(pheromones_0(i+1,j)-pheromones_0(i-1,j))/(2*dx);
    Fy_1(i,j) = alpha1*(pheromones_0(i,j+1)-pheromones_0(i,j-1))/(2*dy);
    
    Fx_1(Fx_1>maxF1)=maxF1;
    Fy_1(Fy_1>maxF1)=maxF1;
    Fx_1(Fx_1<-maxF1)=-maxF1;
    Fy_1(Fy_1<-maxF1)=-maxF1;
    
    ants_0(i,j)=ants_0n(i,j)+(vis_0*dt*(ants_0n(i+1,j)-2*ants_0n(i,j)+ants_0n(i-1,j))/(dx*dx))+(vis_0*dt*(ants_0n(i,j+1)-2*ants_0n(i,j)+ants_0n(i,j-1))/(dy*dy))...
        - dt*zeta_0^(-1)*(FObsx(i,j).*(ants_0n(i+1,j)-ants_0n(i-1,j))/(2*dx)  +  ants_0n(i,j).*(FObsx(i+1,j)-FObsx(i-1,j))/(2*dx))...
        - dt*zeta_0^(-1)*(FObsy(i,j).*(ants_0n(i,j+1)-ants_0n(i,j-1))/(2*dy)  +  ants_0n(i,j).*(FObsy(i,j+1)-FObsy(i,j-1))/(2*dy))...
        - dt*zeta_0^(-1)*(Fx_0(i,j).*(ants_0n(i+1,j)-ants_0n(i-1,j))/(2*dx)  +  ants_0n(i,j).*(Fx_0(i+1,j)-Fx_0(i-1,j))/(2*dx))...
        - dt*zeta_0^(-1)*(Fy_0(i,j).*(ants_0n(i,j+1)-ants_0n(i,j-1))/(2*dy)  +  ants_0n(i,j).*(Fy_0(i,j+1)-Fy_0(i,j-1))/(2*dy)) -...
        dt*loading_rate*food_centered_gaussian(i,j).*ants_0n(i,j) + dt*0*nest_centered_gaussian(i,j) + dt*unloading_rate*nest_centered_gaussian(i,j).*ants_1n(i,j);
    
    
    ants_1(i,j)=ants_1n(i,j)+(vis_1*dt*(ants_1n(i+1,j)-2*ants_1n(i,j)+ants_1n(i-1,j))/(dx*dx))+(vis_1*dt*(ants_1n(i,j+1)-2*ants_1n(i,j)+ants_1n(i,j-1))/(dy*dy)) ...
        - dt*zeta_1^(-1)*(FObsx(i,j).*(ants_1n(i+1,j)-ants_1n(i-1,j))/(2*dx)  +  ants_1n(i,j).*(FObsx(i+1,j)-FObsx(i-1,j))/(2*dx))...
        - dt*zeta_1^(-1)*(FObsy(i,j).*(ants_1n(i,j+1)-ants_1n(i,j-1))/(2*dy)  +  ants_1n(i,j).*(FObsy(i,j+1)-FObsy(i,j-1))/(2*dy))...
        - dt*zeta_1^(-1)*(Fx_1(i,j).*(ants_1n(i+1,j)-ants_1n(i-1,j))/(2*dx)  +  ants_1n(i,j).*(Fx_1(i+1,j)-Fx_1(i-1,j))/(2*dx))...
        - dt*zeta_1^(-1)*(Fy_1(i,j).*(ants_1n(i,j+1)-ants_1n(i,j-1))/(2*dy)  +  ants_1n(i,j).*(Fy_1(i,j+1)-Fy_1(i,j-1))/(2*dy))+...
        dt*loading_rate*food_centered_gaussian(i,j).*ants_0n(i,j) - dt*unloading_rate*nest_centered_gaussian(i,j).*ants_1n(i,j);
    
    %     Boundary conditions
    %     Dirichlet:
    ants_0(1,:)=UW;
    ants_0(nx,:)=UE;
    ants_0(:,1)=US;
    ants_0(:,ny)=UN;
    
    %     ants_0(x-0.8>0,:) = 0;
    
    ants_1(1,:)=UW;
    ants_1(nx,:)=UE;
    ants_1(:,1)=US;
    ants_1(:,ny)=UN;
    
%     ants_1(abs(x-0.6)<0.6,abs(y-1.0)<0.05)=0;
%     ants_0(abs(x-0.6)<0.6,abs(y-1.0)<0.05)=0;
    
    ants_1(U<0.06)=0;
    ants_0(U<0.06)=0;
    %     Neumann:
        ants_0(1,:)=ants_0(2,:)-UnW*dx;
    %     ants_0(nx,:)=ants_0(nx-1,:)+UnE*dx;
    %     ants_0(:,1)=ants_0(:,2)-UnS*dy;
    %     ants_0(:,ny)=ants_0(:,ny-1)+UnN*dy;
    %
    %     ants_1(1,:)=ants_1(2,:)-UnW*dx;
    %     ants_1(nx,:)=ants_1(nx-1,:)+UnE*dx;
    %     ants_1(:,1)=ants_1(:,2)-UnS*dy;
    %     ants_1(:,ny)=ants_1(:,ny-1)+UnN*dy;
    %     }
    if it<nt
        if mod(it,200) == 0
            clf(figure(1));
%             delete(h2);
        end
    end
    % delete(h3);
    
end

if create_movie
    close(mov2);
end