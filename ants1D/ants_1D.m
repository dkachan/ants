close all
clear all
clc
%%
%Specifying Parameters
nx=100;               %Number of steps in space(x)
nt=10000;              %Number of time steps 
dt=0.1;             %Width of each time step
dx=100/(nx-1);         %Width of space step
x=-50:dx:50;            %Range of x (0,2) and specifying the grid points
c1=zeros(nx,1);      %Preallocating c1
c1n=zeros(nx,1);     %Preallocating c1n
c2=zeros(nx,1);      %Preallocating c2
c2n=zeros(nx,1);     %Preallocating c2n
F1 = zeros(nx,1);
F2 = zeros(nx,1);
f1= 200;               %coupling 1
f2= 200;               %coupling 2
CL=1;                %Left Dirichlet B.C
CR=1;                %Right Dirichlet B.C
alpha = 0.8;
MU_1 = -20;
SIGMA_1 = 10;
gaussian1 = normpdf(x,MU_1,SIGMA_1)';

MU_2 = 20;
SIGMA_2 = 10;
gaussian2 = normpdf(x,MU_2,SIGMA_2)';


%% Initial conditions
c1 = gaussian1;
% c2 = gaussian2;

%%
%B.C vector
bc=zeros(nx-2,1);
bc(1)=dt*CL/dx^2; bc(nx-2)=dt*CR/dx^2;  %Dirichlet B.Cs

%%
%Calculating the velocity profile for each time step
i=2:nx-1;
figure
hold on
for it=0:nt
    c1n=c1;
    c2n=c2;
    
    if mod(it,1)==0
        plot(x,c1); 
        hold on
        plot(x,c2);
        plot(x,c2+c1);
        drawnow;
        clf(figure(1))
    end
    
    F1(i) = (c2(i+1)-c2(i-1))/(2*dx);
    F2(i) = (c1(i+1)-c1(i-1))/(2*dx);
    
    c1(i)=c1n(i)+(dt*(c1n(i+1)-2*c1n(i)+c1n(i-1))/(dx*dx)) - dt*f1*(F1(i)-alpha*F2(i)).*(c1n(i+1)-c1n(i-1))/(2*dx) - dt*f1*c1n(i).*(F1(i+1)-F1(i-1)-alpha*(F2(i+1)-F2(i-1)))/(2*dx)...
        + dt*gaussian1(i).*c2(i) - dt*gaussian2(i).*c1(i) ;
    c2(i)=c2n(i)+(dt*(c2n(i+1)-2*c2n(i)+c2n(i-1))/(dx*dx)) - dt*f2*(F2(i)-alpha*F1(i)).*(c2n(i+1)-c2n(i-1))/(2*dx) - dt*f2*c2n(i).*(F2(i+1)-F2(i-1)-alpha*(F1(i+1)-F1(i-1)))/(2*dx)...
        - dt*gaussian1(i).*c2(i) + dt*gaussian2(i).*c1(i) ;
    
    c1(1) = c1(2);
    c1(nx) = c1(nx-1);
    c2(1) = c2(2);
    c2(nx) = c2(nx-1);
  
   
end