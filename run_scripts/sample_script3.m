clc
p = SimParams();
p.N_time_steps = 400000;
p.noise_strength = 15;
p.show_figure = true;

p.nest_center =  [0.25 0.25];
p.init_ant_variance = 0.1;
p.nest_radius = 0.06;
p.ant_velocity = 0.0001;
p.max_pheromone = 20;
p.ant_number = 1000;
p.display_interval =20;
p.beta = 0.002;
p.trail_angle = 0;
sim_ants2(p);
