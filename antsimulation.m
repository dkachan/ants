% ********************************************************************
% *****                     Ant simulation                       *****
% ********************************************************************

% Simple code for ant simulation
% Element Method, FEM

% ---------------------------------------
clear all
clc
tic;
close all
% ---------------------------------------

% ++++++++++++++++++++++
%   GLOBAL VARIABLES
% ++++++++++++++++++++++

global Lx Ly element_length node element topNodes botNodes rightNodes leftNodes numnode center_node dx Edy


% +++++++++++++++++++++++++++++
%             INPUT
% +++++++++++++++++++++++++++++

N_element_y = 30;
N_element_x = 35;

element_length = 1;

Ly = N_element_y*element_length/2 ; %y
Lx = N_element_x*element_length/2; %x

N_time_steps = 10;




% +++++++++++++++++++++++++++++
%       CONSTRUCT MESH
% +++++++++++++++++++++++++++++

% Number of nodes along two directions
nnx = N_element_x+1;
nny = N_element_y+1;

% Four corner points
pt1 = [-Lx/2 -Ly/2] ;pt2 = [Lx/2 -Ly/2] ;pt3 = [Lx/2 Ly/2] ;pt4 = [-Lx/2 Ly/2] ;

% 4-nodes element mesh for F & J
elemType = 'Q4' ;   %???What does this do?
[node,element,~] = meshRectangularRegion(pt1, pt2, pt3, pt4, nnx,nny,elemType);

dx = node(2,1)- node(1,1);
dy = node(1+nnx,2)- node(1,2);

% compute number of nodes, of elements
numnode = size(node,1);
numelem = size(element,1);

% +++++++++++++++++++++++++++++
%       BOUNDARY NODES
% +++++++++++++++++++++++++++++

kt=0;
kb=0;
kl=0;
kr=0;
for i=1:numnode
    if abs(node(i,2)-Ly/2)<1e-10
        kt=kt+1;
        topNodes(kt)=i;
    end
    if abs(node(i,2)+Ly/2)<1e-10;
        kb=kb+1;
        botNodes(kb)=i;
    end
    if abs(node(i,1)+Lx/2)<1e-10;
        kl=kl+1;
        leftNodes(kl)=i;
    end
    if abs(node(i,1)-Lx/2)<1e-10;
        kr=kr+1;
        rightNodes(kr)=i;
    end
end

total_unknown = numnode;

W=zeros(numelem,4);
Q=zeros(numelem,4,2);
center_node = zeros(numelem,2);
for iel=1:numelem
    %find the gauss points coord and weight in 2D regular parent element9
    center_node(iel,:) = mean(node(element(iel,:),:));
    order = 2 ;
    [W(iel,1:order^2),Q(iel,1:order^2,:)] = quadrature(order,'GAUSS',2);
    ngp = 4;
end

% +++++++++++++++++++++++++++++
%             ANT INITIAL CONDITION
% +++++++++++++++++++++++++++++

%inititialize ants position
ant_number = 100;
ant_nest = [-5 0];
ant_pos_center = [ant_nest(1)*ones(ant_number,1),ant_nest(2)*ones(ant_number,1)];
ant_pos_box_size = 2;
init_ant_pos(1:ant_number,1:2) = ant_pos_center+(rand(ant_number,2)-0.5)*ant_pos_box_size;

%inititialize pheronomes
pheronomes = zeros(numnode,ant_number);

%initialize food
init_food_pos = [4 0];
init_food_radius = 2;






figure
hold on
plot_mesh(node,element,elemType,'b-');

figure(2)
figure(1)
%time loop

ant_pos(1:ant_number,1:2) = init_ant_pos;
active_ant = zeros(ant_number,1);
active_trail = zeros(ant_number,1);
dt= 0.8;
ant_velocity=0.15;
evap_speed = 0.001;
nplot =0;
pheronome_gradient = zeros(numelem,2);
active_elements  = zeros(numelem,1);

for it = 1:N_time_steps
    
    nplot = nplot +1;
    
    distance_ant_food =((ant_pos(:,1)-init_food_pos(1)).^2+(ant_pos(:,2)-init_food_pos(2)).^2).^(1/2);
    
    for ia = 1:ant_number
        
        d =((ant_pos(ia,1)-center_node(:,1)).^2+(ant_pos(ia,2)-center_node(:,2)).^2).^(1/2);
        [distances, IDX] = sort(d);
        close_elements = IDX(find(distances<0.5*element_length));
        
        
        iel = IDX(1);
        
        sctr = element(iel,:);
        [N,dNdxi]=lagrange_basis(elemType,[0 0]);
        J0 = node(sctr,:)'*dNdxi;   % element Jacobian matrix
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;
        
        pheron_old = 0;
        for i = 1 : length(active_trail)
            if active_trail(i)==1
                pheron = pheronomes(sctr,i)'*N;
                pheron_gradient = -15*pheron*pheronomes(sctr,i)'*dNdx;
                if pheronomes(sctr(1),i)*pheronomes(sctr(2),i)*pheronomes(sctr(3),i)*pheronomes(sctr(4),i)==0
                    pheron_gradient = [0 0];
                end
                if pheron>pheron_old
                pheronome_gradient(iel,:) = pheron_gradient;
                end
                pheron_old = pheron;
            end
        end
        
        distance_ant_pheronome =((ant_pos(ia,1)-node(unique(element(close_elements,:)),1)).^2+(ant_pos(ia,2)-node(unique(element(close_elements,:)),2)).^2).^(1/2);
        ant_direction(ia,1:2) = (rand(1,2)-0.5);
        
        if norm(ant_pos(ia,1:2)-ant_nest)<0.5
            active_ant(ia) = 0;
        end
        
        
        if distance_ant_food(ia)<init_food_radius || active_ant(ia) ~= 0;
            if active_ant(ia)==0
                active_ant(ia) = ia;
                active_trail(ia) = 1;
                init_food_radius = init_food_radius;
            end
            pheronomes(unique(element(close_elements,:)),ia)=0.5*ones(length(unique(element(close_elements,:))),1);
            ant_direction(ia,1:2) = ant_direction(ia,1:2) - [ant_pos(ia,1:2)-ant_nest]/norm([ant_pos(ia,1:2)-ant_nest]);
        else
            ant_direction(ia,1:2) = ant_direction(ia,1:2) + pheronome_gradient(iel,:)*10;
        end
        
        if abs(ant_pos(ia,1)-(-Lx/2))<1
           ant_direction(ia,1) = 0.5;
        end
        if abs(ant_pos(ia,2)-(Ly/2))<1
           ant_direction(ia,2) = -0.5;
        end
        if abs(ant_pos(ia,1)-(Lx/2))<1
           ant_direction(ia,1) = -0.5;
        end
        if abs(ant_pos(ia,2)-(-Ly/2))<1
           ant_direction(ia,2) = +0.5;
        end
        
        ant_direction(ia,1:2) = ant_direction(ia,1:2)/norm(ant_direction(ia,1:2));
    end
    
    pheronomes = pheronomes - pheronomes.*ones(numnode,ant_number)*evap_speed*dt;
    
    ant_pos(1:ant_number,1:2) = ant_pos(1:ant_number,1:2) + ant_velocity*ant_direction(1:ant_number,1:2)*dt;
    
    h1 = plot(ant_pos(:,1),ant_pos(:,2),'ro');
    
    if any(active_ant)
    h2 = plot(ant_pos(active_ant(active_ant~=0),1),ant_pos(active_ant(active_ant~=0),2),'gx','LineWidth',3);
    end
    h3 = circle(init_food_pos(1),init_food_pos(2),init_food_radius);
    
    h4 = quiver(center_node(:,1),center_node(:,2),pheronome_gradient(:,1),pheronome_gradient(:,2),'black');
    
    pause(0.05)
    
   
    if it<N_time_steps
    delete(h1)
    if any(active_ant)
    delete(h2)
    end
    delete(h3)
    delete(h4)
    end
    
    if it ==1 || nplot == 20
       nplot =0;
       clf(figure(2))
       figure(2)
       hold on
       plot_field(node,element,elemType,sum(pheronomes,2))
       plot_mesh(node,element,elemType,'b-');
       plot(ant_pos(:,1),ant_pos(:,2),'ro');
       plot(ant_pos(active_ant(active_ant~=0),1),ant_pos(active_ant(active_ant~=0),2),'gx','LineWidth',3);
       circle(init_food_pos(1),init_food_pos(2),init_food_radius);
       quiver(center_node(:,1),center_node(:,2),pheronome_gradient(:,1),pheronome_gradient(:,2),'black');
       caxis([0 max(sum(pheronomes,2))]);
    end
    
    figure(1)
    
end