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

left_elem = 1:N_element_x:N_element_y*N_element_x;
bot_elem = 1:N_element_x;
right_elem = N_element_x:N_element_x:N_element_y*N_element_x;
top_elem = N_element_y*(N_element_x-1)+1:N_element_y*N_element_x;

corner1 = left_elem(1);
corner2 = bot_elem(end);
corner3 = right_elem(end);
corner4 = top_elem(1);

left_elem = left_elem(2:end);
right_elem = right_elem(2:end);
top_elem = top_elem(2:end);
bot_elem = bot_elem(2:end);

%find the 8 neighbouring elements [right-middle, right-top,
    % middle-top, left-top, left-middle, left-bottom, middle-bottom,
    % right-bottom]
    
element_list = (1:numelem)';
element_neighbours = [element_list+1, element_list+1+N_element_x, element_list+N_element_x, ...
    element_list+N_element_x-1, element_list-1, element_list-N_element_x-1, element_list-N_element_x,...
    element_list-N_element_x+1];
% element_neighbours(element_neighbours<0) = element_neighbours(element_neighbours<0) + N_element_x;
% Adjust for periodic boundary conditions


if p.periodic_boundaries 
    element_neighbours(bot_elem, 6:8) = [element_list(bot_elem)-1 element_list(bot_elem) element_list(bot_elem)+1] + (N_element_x-1)*N_element_y;
    element_neighbours(right_elem, [1 2 8]) = [element_list(right_elem) element_list(right_elem)+N_element_x element_list(right_elem)-N_element_x]-(N_element_x-1);
    element_neighbours(top_elem, 2:4) = [element_list(top_elem)+1 element_list(top_elem) element_list(top_elem)-1] - (N_element_x-1)*N_element_y;
    element_neighbours(left_elem, 4:6) = [element_list(left_elem)+N_element_x element_list(left_elem) element_list(left_elem)-N_element_x]+(N_element_x-1);
    
    element_neighbours(corner1,:) = [2, N_element_x+2, N_element_x+1, 2*N_element_x,...
        N_element_x, N_element_y*N_element_x, (N_element_y-1)*N_element_x+1, (N_element_y-1)*N_element_x+2];
    element_neighbours(corner2,:) = [1 N_element_x+1, 2*N_element_x, 2*N_element_x-1,...
        N_element_x-1,N_element_y*N_element_x-1, N_element_y*N_element_x, (N_element_y-1)*N_element_x+1];
    element_neighbours(corner3,:) = [(N_element_y-1)*N_element_x+1, 1, N_element_x,...
        N_element_x-1, N_element_y*N_element_x-1,(N_element_y-1)*N_element_x-1,(N_element_y-1)*N_element_x,(N_element_y-2)*N_element_x+1];
    element_neighbours(corner4,:) = [(N_element_y-1)*N_element_x+2, 2, 1,...
        N_element_x, N_element_y*N_element_x, (N_element_y-1)*N_element_x, (N_element_y-2)*N_element_x+1,(N_element_y-2)*N_element_x+2];
    
end

% 
% element_neighbours(1,:)
% element_neighbours(100,:)
% element_neighbours(10000,:)
% element_neighbours(9901,:)
% pause
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
food_boundary_radius = p.food_boundary_radius;
food_boundary_center = p.food_boundary_center;
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
%beta = p.beta;
%nest_drift_velocity = p.nest_drift_velocity;
noise_strength = p.noise_strength;
%obstacle_strength = p.obstacle_strength;
% +++++++++++++++++++++++++++++
%        PHEROMONE MODEL
% +++++++++++++++++++++++++++++
evap_rate = p.evap_rate;
deposition_rate = p.deposition_rate;
release_delay = p.release_delay;


% +++++++++++++++++++++++++++++
%    PREALLOCATION
% +++++++++++++++++++++++++++++

pheromones = zeros(numnode,1);
ant_pos = init_ant_pos;
delayed_ant_pos =zeros(ant_number,2*release_delay);
ant_has_food = false(ant_number,1);
ant_found_food = false(ant_number,1);
prev_angle = 2*pi*rand(ant_number, 1);
delta_theta = zeros(ant_number,1);
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
    clf(1)
    hold on
    if show_mesh
        plot_mesh_dk(node,element,'b-');
    end
    if food_radius
        circle(food_center(1),food_center(2),food_radius);
    end
    if food_boundary_radius
        circle(food_boundary_center(1),food_boundary_center(2),food_boundary_radius);
    end
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
    xlim([0,Lx]);
    ylim([0,Ly]);
    %axis equal
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
    
    if it > release_delay
        delayed_ant_elements_num = floor(delayed_ant_pos(:,1)/element_length) + ... 
        1 + floor(delayed_ant_pos(:,2)/element_length) * N_element_x;
    end
    
    %calculate pheronome levels each ant sees. 
    if ant_number == 1
        ant_pheromones = pheromones(ant_elements)'*N;
        pheromone_gradient = pheromones(ant_elements)'*dNdx;
    else
        ant_pheromones = pheromones(ant_elements)*N;
        pheromone_gradient = pheromones(ant_elements)*dNdx;
    end
    
    %Feed the ants
    if food_radius
        distance_ant_food = sqrt(sum((ant_pos - ones(ant_number,1)*food_center).^2, 2)); 
        near_food =  (distance_ant_food <= food_radius);
        ant_found_food = ~ant_found_food & near_food & ~ant_has_food;
        ant_has_food = ant_has_food | near_food;
    end
    if food_boundary_radius
        distance_ant_food = sqrt(sum((ant_pos - ones(ant_number,1)*food_boundary_center).^2, 2)); 
        near_food = (distance_ant_food >= food_boundary_radius);
        ant_found_food = ~ant_found_food & near_food & ~ant_has_food;
        ant_has_food = ant_has_food | near_food ;
    end
        

    %Drop off food at the nest
    vec_ant_nest = ant_pos - ones(ant_number,1)*nest_center;
    distance_ant_nest = sqrt(sum((vec_ant_nest).^2, 2)); 
    ant_dropped_food = ant_has_food & (distance_ant_nest <= nest_radius);
    ant_has_food = ant_has_food & ~(distance_ant_nest <= nest_radius);
    
    %Dynamics.
    %----angular noise----
    angle_max = p.angle_max;
    tanh_factor = tanh(p.beta*ant_pheromones);
    delta_theta(:) = 0;
    noise = noise_strength*(rand(ant_number, 1) - .5);
    %delta_theta(~ant_has_food) = noise(~ant_has_food);
    delta_theta = noise_strength*(1-tanh_factor).*(rand(ant_number, 1) - .5);
    %----gradient coupling----
    %take 3
    prev_direction = [cos(prev_angle) sin(prev_angle)];
    
    
    n_pheromone_gradient = normalize(pheromone_gradient);
    n_pheromone_gradient(isnan(n_pheromone_gradient)) = 0;
    cp = prev_direction(:,1).*n_pheromone_gradient(:,2) - prev_direction(:,2).*n_pheromone_gradient(:,1);
    %cp = sign(cp).*min(abs(cp),1);
    dot_product = sum(prev_direction.*n_pheromone_gradient,2);
    %dot_product = sign(dot_product).*min(abs(dot_product),1);
    %disp(cp)
    delta_theta = delta_theta + gradient_coupling*(.5*cp - p.gamma*dot_product.*cp);
    
    
    
    %take 2
    %prev_direction = [cos(prev_angle) sin(prev_angle)];
    %dot_product = sum(prev_direction.*pheromone_gradient);
    %summed = prev_direction + pheromone_gradient;
    %angles = atan2(summed(:,2),summed(:,1));
    %disp(dot_product);
    %mask = dot_product > 0;
    %delta_theta(mask) = delta_theta(mask) + (angles(mask) - prev_angle(mask));
    
    %take 1
    %normalized_gradient= normalize(pheromone_gradient);
    %gradient_angle = atan(normalized_gradient(:,2)./normalized_gradient(:,1));
    %phi = gradient_angle - prev_angle;
    %phi(isnan(phi)) = 0;
    %delta_theta = delta_theta + gradient_coupling*sign(phi).*exp(-(abs(phi)-1.5708).^2);
    
    %gradient_coupling*pheromone_gradient;
    
    if it < N_time_steps
        delta_theta = sign(delta_theta).*min(abs(delta_theta),angle_max);
        new_angle = prev_angle + delta_theta;
        %----Turn around after finding food------
        if p.reverse_at_food
            new_angle(ant_found_food) = prev_angle(ant_found_food) + pi;
        end
        if p.reversal_after_food
            new_angle(ant_dropped_food) = prev_angle(ant_dropped_food) + pi;
        end
        if it <= release_delay
            delayed_ant_pos(:,2*it - 1:2*it) = ant_pos;
        else
            delayed_ant_pos = circshift(delayed_ant_pos,[0 -2]);
            delayed_ant_pos(:,end-1:end) = ant_pos;
        end
        ant_pos = ant_pos + ant_velocity*[cos(new_angle) sin(new_angle)]*dt;
        prev_angle = mod(new_angle,2*pi);
    end
    
    
    
    %------------OLD----------------
    %-------FOR REFERENCE-----------
    %TODO ADD OBSTACLES IN NEW FRAMEWORK
    
    %{
    
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
    
   %Obstacle forces

    %
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
    
    
   
    %----Normalize the forces to get a direction-----
    new_ant_direction = normalize(new_ant_direction);
    if it < N_time_steps
        ant_pos = ant_pos + ant_velocity*new_ant_direction*dt;
        prev_direction = new_ant_direction;
    end
    
    %}
    
    %------periodic_bc------
    if p.periodic_boundaries
        ant_pos((ant_pos(:,1)>Lx),1) = ant_pos((ant_pos(:,1)>Lx),1) - Lx;
        ant_pos((ant_pos(:,1)<0),1) = ant_pos((ant_pos(:,1)<0),1) + Lx;
        ant_pos((ant_pos(:,2)>Ly),2) = ant_pos((ant_pos(:,2)>Ly),2) - Ly;
        ant_pos((ant_pos(:,2)<0),2) = ant_pos((ant_pos(:,2)<0),2) + Ly;
    end
    %------Activate pheremone production.  Simply on if you have food, else off-----
    if it>release_delay
    pheromone_elements = delayed_ant_elements_num(ant_has_food); 
    pheromone_elements_withNeighbours = [delayed_ant_elements_num(ant_has_food) element_neighbours(delayed_ant_elements_num(ant_has_food),:)]; 
    
   
    
    if any(pheromone_elements_withNeighbours)
        for j = 1:9
            for i=1:4
                r = sqrt((delayed_ant_pos(ant_has_food,1)-node(element(pheromone_elements_withNeighbours(:,j),i),1)).^2 ...
        + (delayed_ant_pos(ant_has_food,2)-node(element(pheromone_elements_withNeighbours(:,j),i),2)).^2);
                pheromones(element(pheromone_elements,i)) =  pheromones(element(pheromone_elements,i)) + deposition_rate./r;
            end
        end
    end
    
%     if any(pheromone_elements)
%         for i=1:4
%             pheromones(element(pheromone_elements,i)) =  pheromones(element(pheromone_elements,i)) + deposition_rate;
%         end
%     end
    %-----maximum pheromone level-----
    if p.max_pheromone
        pheromones(pheromones > p.max_pheromone) = p.max_pheromone;
    end
    
    %-----free production of pheromones without food
    for i=1:4  
        pheromones(element(delayed_ant_elements_num(~ant_has_food),i)) =  pheromones(element(delayed_ant_elements_num(~ant_has_food),i)) + p.free_deposition_rate;
    end
    %Pheromone volatility
    pheromones = pheromones - pheromones.*ones(numnode,1)*evap_rate*dt;
    end
    
    %----------Real time visualization-----------
    if  show_figure && mod(it,display_interval) == 0
        ant_loc_plot = plot(ant_pos(:,1),ant_pos(:,2),'ro');
        xlim([0,Lx]);
    ylim([0,Ly]);
        if any(ant_has_food)
            ants_with_food_plot =plot(ant_pos(ant_has_food,1),ant_pos(ant_has_food,2),'gx','LineWidth',3);
            xlim([0,Lx]);
    ylim([0,Ly]);
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
            xlim([0,Lx]);
    ylim([0,Ly]);
            if any(ant_has_food) 
                delete(ants_with_food_plot)
                xlim([0,Lx]);
    ylim([0,Ly]);
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