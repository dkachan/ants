% ********************************************************************
% *****                     Ant simulation                       *****
% ********************************************************************

% Simple code for ant simulation
% Element Method, FEM

% ---------------------------------------
%clear all
clc
tic;
%close all
% ---------------------------------------

% ++++++++++++++++++++++
%   GLOBAL VARIABLES
% ++++++++++++++++++++++

global Lx Ly element_length node element numnode center_node 


% +++++++++++++++++++++++++++++
%             INPUT
% +++++++++++++++++++++++++++++

N_element_y = 100;
N_element_x = 100;
element_length = 1;

Ly = N_element_y*element_length ; %y
Lx = N_element_x*element_length; %x

N_time_steps = 2000;


% +++++++++++++++++++++++++++++
%       CONSTRUCT MESH
% +++++++++++++++++++++++++++++

% Number of nodes along two directions
nnx = N_element_x+1;
nny = N_element_y+1;

% Four corner points
pt1 = [0 0] ;pt2 = [Lx 0] ;pt3 = [Lx Ly] ;pt4 = [0 Ly] ;

% 4-nodes element mesh for F & J
elemType = 'Q4' ;   
[node,element,~] = meshRectangularRegion(pt1, pt2, pt3, pt4, nnx,nny,elemType);

numnode = size(node,1);
numelem = size(element,1);

center_node = zeros(numelem,2);
for iel=1:numelem
    center_node(iel,:) = mean(node(element(iel,:),:));
end


nFrames = 250;
frame_interval = floor(N_time_steps/nFrames);
mov(1:nFrames) = struct('cdata',[], 'colormap',[]);
k=1;




% +++++++++++++++++++++++++++++
%             ANT INITIAL CONDITION
% +++++++++++++++++++++++++++++
averages = zeros(20,1);
distances = linspace(70,80,5);
for a = 1:1
    for b=1:1
%inititialize ants position
ant_number = 6000;
ant_nest = [20 50];
ant_pos_center = [ant_nest(1)*ones(ant_number,1),ant_nest(2)*ones(ant_number,1)];
ant_pos_box_size = 2;
init_ant_pos = ant_pos_center+(rand(ant_number,2)-0.5)*ant_pos_box_size;
%inititialize pheronomes
pheromones_dk = zeros(numnode,1);

%initialize food
init_food_pos = [distances(a) 50];  %73 is critical!
init_food_radius = 5;
init_nest_radius = 5;

show_first_fig = true;
show_second_fig = false;

obstacle_center = [48 8];
obstacle_radius = 8;


if show_first_fig
    figure(1)
    hold on
    %plot_mesh_dk(node,element,'b-');
    circle(init_food_pos(1),init_food_pos(2),init_food_radius);
    circle(ant_nest(1),ant_nest(2),init_nest_radius);
    ho = circle(obstacle_center(1),obstacle_center(2),obstacle_radius);
    axis([0,100,0,100])
end
if show_second_fig
    figure(2)
    plot_mesh_dk(node,element,'b-');
    circle(init_food_pos(1),init_food_pos(2),init_food_radius);
    circle(ant_nest(1),ant_nest(2),init_nest_radius);
    figure(1)
end

ant_pos(1:ant_number,1:2) = init_ant_pos;
ant_has_food = false(ant_number,1);
ant_releasing_pheromones = false(ant_number, 1);
prev_direction = 1/sqrt(2) * ones(ant_number, 2);
dt= 1;
ant_velocity=1.075;
evap_speed = 0.001;
deposition_rate = .02;
alpha = .15;
beta = 1;
drift_velocity = .4;



[N,dNdxi]=lagrange_basis(elemType,[0 0]);
dNdx  = dNdxi*inv(node(element(1,:),:)'*dNdxi);


food_trace = zeros(N_time_steps,1);
for it = 1:N_time_steps
        
    distance_ant_food = ((ant_pos(:,1)-init_food_pos(1)).^2+(ant_pos(:,2)-init_food_pos(2)).^2).^(1/2);
    distance_ant_nest = ((ant_pos(:,1)-ant_nest(1)).^2+(ant_pos(:,2)-ant_nest(2)).^2).^(1/2);
    distance_ant_obstacle = ((ant_pos(:,1)-obstacle_center(1)).^2+(ant_pos(:,2)-obstacle_center(2)).^2).^(1/2);

    %find elements of ants
    ant_elements = element(floor(ant_pos(:,1)/element_length) + 1 + floor(ant_pos(:,2)/element_length) * N_element_x, :);
    
    %calculate pheronome terms, the pheromone level that each ant sees

    pheromone = pheromones_dk(ant_elements)*N;
    pheromone_gradient = [pheromones_dk(ant_elements)*dNdx(:,1) pheromones_dk(ant_elements)*dNdx(:,2)] ;
    
    %Feed the ants
    ant_has_food = ant_has_food | (distance_ant_food <= init_food_radius);
    
    %Drop off food at the nest
    dropped_food = ant_has_food & (distance_ant_nest <= init_nest_radius);
    ant_has_food = ant_has_food & ~(distance_ant_nest <= init_nest_radius);
    
    %Dynamics
    noise = rand(ant_number,2)-0.5;
    drift = -drift_velocity*(ant_has_food*[1 1]).*[(ant_pos(:,1) - ant_nest(1)) (ant_pos(:,2) - ant_nest(2))]./(distance_ant_food*[1 1]);
    pheromone_force = alpha*(~ant_has_food*[1 1]).*pheromone_gradient;
    
    near_obstacle = (distance_ant_obstacle <= obstacle_radius);
    obstacle_force = 5*(near_obstacle*[1 1]).*[(ant_pos(:,1) - obstacle_center(1)) (ant_pos(:,2) - obstacle_center(2))]./(distance_ant_obstacle*[1 1]);
    
    
    nest_force = -5*(dropped_food*[1 1]).*prev_direction;
    
    ant_direction = normalize((1-tanh(beta*pheromone)*[1 1]).*noise+nest_force+obstacle_force+drift+pheromone_force + (tanh(beta*pheromone)*[1 1]).*prev_direction);
    if it ~= N_time_steps
        ant_pos = ant_pos + ant_velocity*ant_direction*dt;
        prev_direction = ant_direction;
    end
    
    
    %update obstacle
    obstacle_center(2) = obstacle_center(2) + 80/N_time_steps;
    
    %periodic_bc
    ant_pos((ant_pos(:,1)>Lx),1) = ant_pos((ant_pos(:,1)>Lx),1) - Lx;
    ant_pos((ant_pos(:,1)<0),1) = ant_pos((ant_pos(:,1)<0),1) + Lx;

    ant_pos((ant_pos(:,2)>Ly),2) = ant_pos((ant_pos(:,2)>Ly),2) - Ly;
    ant_pos((ant_pos(:,2)<0),2) = ant_pos((ant_pos(:,2)<0),2) + Ly;
    %Activate pheremone production.  For now this is simply on if you have
    %food, else off. 
    ant_releasing_pheromones = ant_has_food;
    
    %Pheromone release.  Every node acquires the same amount of pheromone.
    pheromone_elements = ant_elements(ant_releasing_pheromones,:); 
    if any(pheromone_elements)
        for i=1:4
            pheromones_dk(pheromone_elements(:,i)) =  pheromones_dk(pheromone_elements(:,i)) + deposition_rate;
        end
    end

    %Pheromone volatility
    pheromones_dk = pheromones_dk - pheromones_dk.*ones(numnode,1)*evap_speed*dt;
    
    
    %if it == N_time_steps && show_first_fig
    if  show_first_fig
        figure(1)
        h1 = plot(ant_pos(:,1),ant_pos(:,2),'ro');
        delete(ho)
        ho =circle(obstacle_center(1),obstacle_center(2),obstacle_radius);

        if any(ant_has_food)
            %pheromone_gradient_el = [pheromones_dk(element)*dNdx(:,1) pheromones_dk(element)*dNdx(:,2)] ;
            %h2 = quiver(center_node(:,1),center_node(:,2),pheromone_gradient_el(:,1),pheromone_gradient_el(:,2),'black');
            h3 =plot(ant_pos(ant_has_food,1),ant_pos(ant_has_food,2),'gx','LineWidth',3);
        end     
        if mod(it,frame_interval) ==0
           %mov(k) = getframe(gca);
           k=k+1;
        end
        pause(0.004)
        if it<N_time_steps
            delete(h1)
            if any(ant_has_food) 
                %delete(h2)
                delete(h3)
            end
        end
    end
    

    if mod(it,100) == 0 && show_second_fig
       figure(2)
       clf()
       hold on
       plot_field(node,element,elemType,pheromones_dk)
       plot(ant_pos(:,1),ant_pos(:,2),'ro');
       plot(ant_pos(ant_has_food,1),ant_pos(ant_has_food,2),'gx','LineWidth',3);
       circle(init_food_pos(1),init_food_pos(2),init_food_radius);
       pheromone_gradient_el = [pheromones_dk(element)*dNdx(:,1) pheromones_dk(element)*dNdx(:,2)] ;
       quiv = quiver(center_node(:,1),center_node(:,2),pheromone_gradient_el(:,1),pheromone_gradient_el(:,2),'black');
       caxis([0 max(pheromones_dk)]);
    end
    
    
food_trace(it) = (sum(ant_has_food) / ant_number);   
end

averages(a) = averages(a) + (sum(ant_has_food) / ant_number)
    end
end
%# save as AVI file, and open it using system video player
%movie2avi(mov, 'moving_obstacle_ants.avi', 'compression','None', 'fps',10);
toc