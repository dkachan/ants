p = SimParams();
p.N_time_steps = 4000;
p.display_interval = 100;
p.nest_drift_velocity = 0;
p.noise_strength = 1.4;
p.show_quiver = false;
p.show_figure = true;
%Obstacles!
o1 = LineObstacle([40 40], [40 60], 1, ObstacleTypes.DoubleVerticalLine);
o2 = LineObstacle([50 40], [80 40], 1, ObstacleTypes.DoubleHorizontalLine);
o3 = LineObstacle([50 60], [80 60], 1, ObstacleTypes.DoubleHorizontalLine);
o4 = CircleObstacle([20 25], 7);

p.obstacles = [ LineObstacle([50 20], [50 60], 2, ObstacleTypes.DoubleVerticalLine) ...
    LineObstacle([70 20], [70 80], 2, ObstacleTypes.DoubleVerticalLine)...
    LineObstacle([30 60], [50 60], 2, ObstacleTypes.DoubleHorizontalLine) ...
    LineObstacle([30 80], [70 80], 2, ObstacleTypes.DoubleHorizontalLine)...
    LineObstacle([50 20], [70 20], 2, ObstacleTypes.DoubleHorizontalLine)];


%p.obstacles = [];

p.obstacle_strength = 5;
p.nest_center =  [30 30];
p.food_center = [60,45];
p.food_radius = 2;
p.nest_radius = 6;
p.ant_velocity = 2.;
p.max_pheromone = 10;
p.evap_rate = .00001;
p.deposition_rate = .0011;
p.gradient_coupling=5.5;
p.ant_number = 3000;
p.free_deposition_rate = 0;
p.beta = .2;
sim_ants(p);
