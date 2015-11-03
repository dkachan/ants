close all
clear
clc

addpath('../aux');
%% Specifying parameters
nx=60;                           %Number of steps in space(x)
ny=60;                           %Number of steps in space(y)
nt=30000;                           %Number of time steps
dt=0.10;                         %Width of each time step
dx=200/(nx-1);                     %Width of space step(x)
dy=200/(ny-1);                     %Width of space step(y)
x=-100:dx:100;                        %Range of x(0,2) and specifying the grid points
y=-100:dy:100;                        %Range of y(0,2) and specifying the grid points
c1=zeros(nx,ny);                  %Preallocating u
c2=zeros(nx,ny);                  %Preallocating u
F1x=zeros(nx,ny);                  %Preallocating Fx_0
F1y=zeros(nx,ny);                  %Preallocating Fy_0
F2x=zeros(nx,ny);                  %Preallocating Fx
F2y=zeros(nx,ny);                  %Preallocating Fy_0
f1 = 800;                      %mobility
f2 = 800;                       %mobility
UW=0;                            %x=0 Dirichlet B.C
UE=0;                            %x=L Dirichlet B.C
US=0;                            %y=0 Dirichlet B.C
UN=0;                            %y=L Dirichlet B.C
UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/dn=UnN)
x1 = -25; %y
y1 = 20; %x
x2 = 25; %y
y2 = -20; %x
alpha = 1;
D1 = 0.2;
D2 = 0.2;
rate1 = 200;
rate2 = 200;

%% Initial Conditions

[X, Y] = meshgrid(x,y);
mu1 = [x1 y1];
mu2 = [x2 y2];
Sigma1 = [10 0; 0 10];
Sigma2 = [10 0; 0 10];
G1 = mvnpdf([X(:) Y(:)],mu1,Sigma1);
G2 = mvnpdf([X(:) Y(:)],mu2,Sigma2);
Gaussian1 = reshape(G1,length(x),length(y));
Gaussian2 = reshape(G2,length(x),length(y));

c1 = reshape(G1,length(x),length(y));
c2 = 0*reshape(G2,length(x),length(y));

create_movie = false;

if create_movie
    mov2 = VideoWriter('ants_smoluchowski4.avi');
    open(mov2);
end

%% Calculating the field variable for each time step



i=2:nx-1;
j=2:ny-1;

hFig = figure('Color',[1 1 1]);
set(hFig, 'Position', [100 50 1000 400])
hold on
plot3(x1,y1,1)
plot3(x2,y2,1)


for it=0:nt
    c1n=c1;
    c2n=c2;
    
    if mod(it,200) == 0
        subplot(1,3,1)
        hold on
        h=surf(x,y,c1,'EdgeColor','none');
        drawnow;
%         view(0,0)
        pause(0.01)
        subplot(1,3,2)
        hold on
        h2=surf(x,y,c2,'EdgeColor','none');
        drawnow;
%         view(0,0)
        pause(0.01)
        subplot(1,3,3)
        hold on
        h2=surf(x,y,c2+c1,'EdgeColor','none');
        drawnow;
%         view(0,0)
        pause(0.01)
        if create_movie
            currFrame = getframe(hFig);
            writeVideo(mov2,currFrame);
        end
    end
    beta =1;
    F1x(i,j) = (c2(i+1,j).^(beta)-c2(i-1,j).^(beta))/(2*dx);
    F1y(i,j) = (c2(i,j+1).^(beta)-c2(i,j-1).^(beta))/(2*dy);
   
    F2x(i,j) = (c1(i+1,j).^(beta)-c1(i-1,j).^(beta))/(2*dx);
    F2y(i,j) = (c1(i,j+1).^(beta)-c1(i,j-1).^(beta))/(2*dy);
    
    c1(i,j)=c1n(i,j)+D1*(dt*(c1n(i+1,j)-2*c1n(i,j)+c1n(i-1,j))/(dx*dx))+D1*(dt*(c1n(i,j+1)-2*c1n(i,j)+c1n(i,j-1))/(dy*dy))...
        - dt*f1*((F1x(i,j)-alpha*F2x(i,j)).*(c1n(i+1,j)-c1n(i-1,j))/(2*dx)  +  c1n(i,j).*(F1x(i+1,j)-F1x(i-1,j)-alpha*(F2x(i+1,j)-F2x(i-1,j)))/(2*dx))...
        - dt*f1*((F1y(i,j)-alpha*F2y(i,j)).*(c1n(i,j+1)-c1n(i,j-1))/(2*dy)  +  c1n(i,j).*(F1y(i,j+1)-F1y(i,j-1)-alpha*(F2y(i,j+1)-F2y(i,j-1)))/(2*dy)) -...
        dt*rate1*Gaussian2(i,j).*c1n(i,j) + dt*rate2*Gaussian1(i,j).*c2n(i,j);
    
    
    c2(i,j)=c2n(i,j)+D2*(dt*(c2n(i+1,j)-2*c2n(i,j)+c2n(i-1,j))/(dx*dx))+D2*(dt*(c2n(i,j+1)-2*c2n(i,j)+c2n(i,j-1))/(dy*dy)) ...
        - dt*f2*((F2x(i,j)-alpha*F1x(i,j)).*(c2n(i+1,j)-c2n(i-1,j))/(2*dx)  +  c2n(i,j).*(F2x(i+1,j)-F2x(i-1,j)-alpha*(F1x(i+1,j)-F1x(i-1,j)))/(2*dx))...
        - dt*f2*((F2y(i,j)-alpha*F1y(i,j)).*(c2n(i,j+1)-c2n(i,j-1))/(2*dy)  +  c2n(i,j).*(F2y(i,j+1)-F2y(i,j-1)-alpha*(F1y(i,j+1)-F1y(i,j-1)))/(2*dy))+...
        dt*rate1*Gaussian2(i,j).*c1n(i,j) - dt*rate2*Gaussian1(i,j).*c2n(i,j);
    
    
    %     Neumann:
    c1(1,:)=c1(2,:)-UnW*dx;
    c1(nx,:)=c1(nx-1,:)+UnE*dx;
    c1(:,1)=c1(:,2)-UnS*dy;
    c1(:,ny)=c1(:,ny-1)+UnN*dy;
    
    c2(1,:)=c2(2,:)-UnW*dx;
    c2(nx,:)=c2(nx-1,:)+UnE*dx;
    c2(:,1)=c2(:,2)-UnS*dy;
    c2(:,ny)=c2(:,ny-1)+UnN*dy;
    
  
    if it<nt
        if mod(it,200) == 0
            clf(figure(1));
%             delete(h2);
        end
    end
    
end

if create_movie
    close(mov2);
end