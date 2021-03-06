
close all
clc

p = SimParams();
p.N_time_steps = 50000;
p.display_interval = 100;
p.ant_number = 200;
p.show_quiver = true;
p.free_deposition_rate = .5;
p.deposition_rate = 2;
p.max_pheromone = 1000;
p.evap_rate = .0005;
p.ant_velocity = .2;
p.pause_time = .002;
p.show_mesh = false;
p.show_figure = true;
p.Lx = 100;
p.Ly = 100;
p.N_element_x=100;
p.N_element_y=100;
p.dt=1;
p.food_boundary_center = [50 50];
p.food_boundary_radius = 0;
p.nest_center = [20 20];
p.food_radius = 10;
p.food_center = [70 70];
p.nest_radius = 10;
p.noise_strength = 1;
p.angle_max = pi/80;
p.reversal_after_food=true;
p.gradient_coupling = .2;
p.release_delay = 15;
p.beta = 100 ;
p.gamma = 8;

sim_ants(p);



%{

This works with a tanh factor on the noise only.  

p = SimParams();
p.N_time_steps = 50000;
p.display_interval = 10;
p.ant_number = 1;
p.show_quiver = true;
p.free_deposition_rate = 1;
p.deposition_rate = 1;
p.max_pheromone = 1000;
p.evap_rate = .01;
p.ant_velocity = 1;
p.pause_time = .002;
p.show_mesh = false;
p.show_figure = true;
p.Lx = 100;
p.Ly = 100;
p.N_element_x=100;
p.N_element_y=100;
p.dt=.2;
p.food_boundary_center = [50 50];
p.food_boundary_radius = 30;
p.nest_center = [50 50];
p.food_radius = 0;
p.food_center = [70 50];
p.nest_radius = 10;
p.noise_strength = 1;
p.angle_max = pi/50;
p.reversal_after_food=true;
p.gradient_coupling = .2;
p.release_delay = 15;
p.beta = 10 ;
p.gamma = 2.5;

%}

%{

great one!!

p = SimParams();
p.N_time_steps = 50000;
p.display_interval = 100;
p.ant_number = 20;
p.show_quiver = true;
p.free_deposition_rate = .5;
p.deposition_rate = 1;
p.max_pheromone = 1000;
p.evap_rate = .0005;
p.ant_velocity = .2;
p.pause_time = .002;
p.show_mesh = false;
p.show_figure = true;
p.Lx = 100;
p.Ly = 100;
p.N_element_x=100;
p.N_element_y=100;
p.dt=1;
p.food_boundary_center = [50 50];
p.food_boundary_radius = 30;
p.nest_center = [50 50];
p.food_radius = 0;
p.food_center = [70 50];
p.nest_radius = 10;
p.noise_strength = 1;
p.angle_max = pi/80;
p.reversal_after_food=true;
p.gradient_coupling = .2;
p.release_delay = 15;
p.beta = 100 ;
p.gamma = 8;


sim_ants(p);
%}

%{

SOLUTION!!!

p = SimParams();
p.N_time_steps = 50000;
p.display_interval = 100;
p.ant_number = 200;
p.show_quiver = true;
p.free_deposition_rate = .5;
p.deposition_rate = 2;
p.max_pheromone = 1000;
p.evap_rate = .0005;
p.ant_velocity = .2;
p.pause_time = .002;
p.show_mesh = false;
p.show_figure = true;
p.Lx = 100;
p.Ly = 100;
p.N_element_x=100;
p.N_element_y=100;
p.dt=1;
p.food_boundary_center = [50 50];
p.food_boundary_radius = 0;
p.nest_center = [20 20];
p.food_radius = 10;
p.food_center = [70 70];
p.nest_radius = 10;
p.noise_strength = 1;
p.angle_max = pi/80;
p.reversal_after_food=true;
p.gradient_coupling = .2;
p.release_delay = 15;
p.beta = 100 ;
p.gamma = 8;

sim_ants(p);


%}