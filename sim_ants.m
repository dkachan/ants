function  out  = sim_ants( p )
%SIM_ANTS Simulates ants with parameters p.
%   Detailed explanation goes here

% ********************************************************************
% *****                     Parameters                           *****
% ********************************************************************

% +++++++++++++++++++++++++++++
%             MESH
% +++++++++++++++++++++++++++++

%------Input--------------
N_element_y = p.N_element_y;
N_element_x = p.N_element_x;
Ly = p.Ly;
Lx = p.Lx; 
elemType = p.elem_type;   

%------Calculate Mesh------
element_length = Lx/N_element_x; %TODO allow for a rectangular mesh.
n_nodes_x = N_element_x+1; 
n_nodes_y = N_element_y+1; 
pt1 = [0 0] ;pt2 = [Lx 0] ;pt3 = [Lx Ly] ;pt4 = [0 Ly] ; %Corner points
[node,element,~] = meshRectangularRegion(pt1, pt2, pt3, pt4, n_nodes_x,n_nodes_y,elemType);
numnode = size(node,1);
numelem = size(element,1);
center_node = [mean(reshape(node(element(:,:),1),numelem,4),2) ...
    mean(reshape(node(element(:,:),2),numelem,4),2) ];

%------Basis vectors & gradient------
[N,dNdxi]=lagrange_basis(elemType,[0 0]);
dNdx  = dNdxi/(node(element(1,:),:)'*dNdxi);

% +++++++++++++++++++++++++++++
%       SIM SETTINGS
% +++++++++++++++++++++++++++++
N_time_steps = p.N_time_steps;
ant_number = p.ant_number;
dt = p.dt;

% +++++++++++++++++++++++++++++
%    GEOMETRY
% +++++++++++++++++++++++++++++
nest_center = p.nest_center;
food_center = p.food_center;  
nest_radius = p.nest_radius;
food_radius = p.food_radius;
%------Obstacles----------
obstacles = p.obstacles;

% +++++++++++++++++++++++++++++
%    INITIAL ANT DISTRIBUTION
% +++++++++++++++++++++++++++++
init_ant_pos = ones(ant_number,1)*nest_center+(rand(ant_number,2)-0.5)*p.init_ant_variance;


% +++++++++++++++++++++++++++++
%            ANT MODEL
% +++++++++++++++++++++++++++++

ant_velocity = p.ant_velocity;
gradient_coupling = p.gradient_coupling;
beta = p.beta;
nest_drift_velocity = p.nest_drift_velocity;
noise_strength = p.noise_strength;
obstacle_strength = p.obstacle_strength;
% +++++++++++++++++++++++++++++
%        PHEROMONE MODEL
% +++++++++++++++++++++++++++++
evap_rate = p.evap_rate;
deposition_rate = p.deposition_rate;


% +++++++++++++++++++++++++++++
%    PREALLOCATION
% +++++++++++++++++++++++++++++

pheromones = zeros(numnode,1);
ant_pos = init_ant_pos;
ant_has_food = false(ant_number,1);
prev_direction = 1/sqrt(2) * ones(ant_number, 2);
%--------output------------
%food_trace = zeros(N_time_steps,1);

% +++++++++++++++++++++++++++++
%    VISUALIZATION
% +++++++++++++++++++++++++++++
show_figure = p.show_figure;
show_mesh = p.show_mesh;
show_quiver = p.show_quiver;
display_interval = p.display_interval;
pause_time = p.pause_time;
create_movie = p.create_movie;


% ********************************************************************
% *****                      MAIN FLOW                           *****
% ********************************************************************

%-----Visualize if necessary---------
if show_figure
    figure(1)
    hold on
    if show_mesh
        plot_mesh_dk(node,element,'b-');
    end
    circle(food_center(1),food_center(2),food_radius);
    circle(nest_center(1),nest_center(2),nest_radius);
    %----Draw obstacles-----
    if ~isempty(obstacles)
        for i = 1:length(obstacles)
            switch obstacles(i).type
                case {ObstacleTypes.VerticalLineLeft, ObstacleTypes.VerticalLineRight,...
                        ObstacleTypes.DoubleVerticalLine, ObstacleTypes.DoubleHorizontalLine, ...
                        ObstacleTypes.HorizontalLineDown, ObstacleTypes.HorizontalLineUp}
                    plot([obstacles(i).pt1(1) obstacles(i).pt2(1)],...
                        [obstacles(i).pt1(2) obstacles(i).pt2(2)],'k','LineWidth',3);
                case ObstacleTypes.Circle
                    circle(obstacles(i).center(1),obstacles(i).center(2),obstacles(i).radius);
            end
        end
    end
    %ho = circle(obstacle_center(1),obstacle_center(2),obstacle_radius);
    axis([0,Lx,0,Ly])
    axis equal
end

if create_movie
    frame_interval = floor(N_time_steps/p.movie_num_frames);
    vidObj = VideoWriter(p.movie_filename);
    vidObj.Quality = p.movie_quality;
    vidObj.FrameRate = p.movie_fps;
    open(vidObj);
    movie_counter=1;
end

%-----Time Iteration--------

for it = 1:N_time_steps
    %find element in which each ant presently resides
    ant_elements = element(floor(ant_pos(:,1)/element_length) + ... 
        1 + floor(ant_pos(:,2)/element_length) * N_element_x, :);
    
    %calculate pheronome levels each ant sees. 
    ant_pheromones = pheromones(ant_elements)*N;
    pheromone_gradient = pheromones(ant_elements)*dNdx;
    
    %Feed the ants
    distance_ant_food = sqrt(sum((ant_pos - ones(ant_number,1)*food_center).^2, 2)); 
    ant_has_food = ant_has_food | (distance_ant_food <= food_radius);

    %Drop off food at the nest
    vec_ant_nest = ant_pos - ones(ant_number,1)*nest_center;
    distance_ant_nest = sqrt(sum((vec_ant_nest).^2, 2)); 
    ant_dropped_food = ant_has_food & (distance_ant_nest <= nest_radius);
    ant_has_food = ant_has_food & ~(distance_ant_nest <= nest_radius);
    
    %Dynamics.
    %----random noise----
    noise = noise_strength*(rand(ant_number,2)-0.5);
    if p.ballistic_coupling
        tanh_factor = tanh(beta*ant_pheromones)*[1 1];
        new_ant_direction = (1-tanh_factor).*noise;
        new_ant_direction = new_ant_direction + tanh_factor.*prev_direction;
    else
        new_ant_direction = noise;
    end
    %----Drift back to the nest with food----
    if p.nest_attraction_with_food
        drift = -nest_drift_velocity*(ant_has_food*[1 1]).*vec_ant_nest./(distance_ant_nest*[1 1]);
        new_ant_direction = new_ant_direction + drift;
    end
    %----Climb up gradients in pheromone concentration----
    if p.gradient_force
        pheromone_force = gradient_coupling*(~ant_has_food*[1 1]).*pheromone_gradient;
        %pheromone_force = gradient_coupling*pheromone_gradient;
        new_ant_direction = new_ant_direction + pheromone_force;
    end 
    %----Turn around after dropping off food------
    if p.reversal_after_food
        nest_force = -p.reversal_strength*(ant_dropped_food*[1 1]).*prev_direction;
        new_ant_direction = new_ant_direction + nest_force;
    end
    
    if ~isempty(obstacles)
        for i = 1:length(obstacles)
            switch obstacles(i).type
                case ObstacleTypes.VerticalLineLeft
                    near_obstacle = (ant_pos(:,2) > obstacles(i).pt1(2)) ...
                        & (ant_pos(:,2) < obstacles(i).pt2(2)) ...
                        & (obstacles(i).pt1(1) - ant_pos(:,1) > 0) ...
                        & (obstacles(i).pt1(1) - ant_pos(:,1) < obstacles(i).width);
                    obstacle_force = obstacle_strength*(near_obstacle*[-1 0]);
                    new_ant_direction = new_ant_direction + obstacle_force;
                        
                case ObstacleTypes.VerticalLineRight
                    near_obstacle = (ant_pos(:,2) > obstacles(i).pt1(2)) ...
                        & (ant_pos(:,2) < obstacles(i).pt2(2)) ...
                        & (ant_pos(:,1) -obstacles(i).pt1(1)  > 0) ...
                        & (ant_pos(:,1) -obstacles(i).pt1(1)  < obstacles(i).width);
                    obstacle_force = obstacle_strength*(near_obstacle*[1 0]);
                    new_ant_direction = new_ant_direction + obstacle_force;
                
                case ObstacleTypes.DoubleVerticalLine
                    near_obstacle = (ant_pos(:,2) > obstacles(i).pt1(2)) ...
                        & (ant_pos(:,2) < obstacles(i).pt2(2)) ...
                        & (abs(ant_pos(:,1) - obstacles(i).pt1(1)) < obstacles(i).width ) ;
                
                    obstacle_force = obstacle_strength*(near_obstacle.*sign(ant_pos(:,1) - obstacles(i).pt1(1)))*[1 0];
                    new_ant_direction = new_ant_direction + obstacle_force;
   
                case ObstacleTypes.DoubleHorizontalLine
                    near_obstacle = (ant_pos(:,1) > obstacles(i).pt1(1)) ...
                        & (ant_pos(:,1) < obstacles(i).pt2(1)) ...
                        & (abs(ant_pos(:,2) - obstacles(i).pt1(2)) < obstacles(i).width ) ;
                
                    obstacle_force = obstacle_strength*(near_obstacle.*sign(ant_pos(:,2) - obstacles(i).pt1(2)))*[0 1];
                    new_ant_direction = new_ant_direction + obstacle_force;
   
                case ObstacleTypes.Circle
                    ant_obstacle_vec = ant_pos - ones(ant_number,1)*obstacles(i).center;
                    distance_ant_obstacle = sqrt(sum(ant_obstacle_vec.^2, 2)); 
                    near_obstacle = (distance_ant_obstacle <= obstacles(i).radius);  
                    obstacle_force = obstacle_strength*(near_obstacle*[1 1]).*ant_obstacle_vec./(distance_ant_obstacle*[1 1]);
                    new_ant_direction = new_ant_direction + obstacle_force;
            end
        end
    end
    
    %Obstacle forces
    
   
    %----Normalize the forces to get a direction-----
    new_ant_direction = normalize(new_ant_direction);
    if it < N_time_steps
        ant_pos = ant_pos + ant_velocity*new_ant_direction*dt;
        prev_direction = new_ant_direction;
    end
    
    %------periodic_bc------
    if p.periodic_boundaries
        ant_pos((ant_pos(:,1)>Lx),1) = ant_pos((ant_pos(:,1)>Lx),1) - Lx;
        ant_pos((ant_pos(:,1)<0),1) = ant_pos((ant_pos(:,1)<0),1) + Lx;
        ant_pos((ant_pos(:,2)>Ly),2) = ant_pos((ant_pos(:,2)>Ly),2) - Ly;
        ant_pos((ant_pos(:,2)<0),2) = ant_pos((ant_pos(:,2)<0),2) + Ly;
    end
    %------Activate pheremone production.  Simply on if you have food, else off-----
    pheromone_elements = ant_elements(ant_has_food,:); 
    if any(pheromone_elements)
        for i=1:4
            pheromones(pheromone_elements(:,i)) =  pheromones(pheromone_elements(:,i)) + deposition_rate;
        end
    end
    %-----maximum pheromone level-----
    if p.max_pheromone
        pheromones(pheromones > p.max_pheromone) = p.max_pheromone;
    end
    
    %-----free production of pheromones
    for i=1:4
        pheromones(ant_elements(:,i)) =  pheromones(ant_elements(:,i)) + p.free_deposition_rate;
    end
    %Pheromone volatility
    pheromones = pheromones - pheromones.*ones(numnode,1)*evap_rate*dt;
    
    %----------Real time visualization-----------
    if  show_figure && mod(it,display_interval) == 0
        ant_loc_plot = plot(ant_pos(:,1),ant_pos(:,2),'ro');
        if any(ant_has_food)
            ants_with_food_plot =plot(ant_pos(ant_has_food,1),ant_pos(ant_has_food,2),'gx','LineWidth',3);
        end     
        if show_quiver
            pheromone_gradient_el = pheromones(element)*dNdx ;
            pheromone_quiver = quiver(center_node(:,1),center_node(:,2),pheromone_gradient_el(:,1),pheromone_gradient_el(:,2),'black');
        end
        if create_movie && mod(it,frame_interval) == 0
           writeVideo(vidObj,getFrame(gca))
           movie_counter=movie_counter+1;
        end
        pause(pause_time)
        if it<N_time_steps
            delete(ant_loc_plot)
            if any(ant_has_food) 
                delete(ants_with_food_plot)
            end
            if show_quiver
                delete(pheromone_quiver)
            end
        end
    end
    
%food_trace(it) = (sum(ant_has_food) / ant_number);   
end


%--------Close video object-------------
if create_movie
    close(vidObj)
end

%placeholder
out = 1;

end