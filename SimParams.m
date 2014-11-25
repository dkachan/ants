classdef SimParams
    %SIM_PARAMS Container for ant simulations parameters
    %   Detailed explanation goes here
    
    properties
        
        % +++++++++++++++++++++++++++++
        %             MESH
        % +++++++++++++++++++++++++++++

        N_element_y = 100;
        N_element_x = 100;
        elem_type = 'Q4';
        
        % +++++++++++++++++++++++++++++
        %           Simulation
        % +++++++++++++++++++++++++++++ 
        
        N_time_steps = 10000;
        dt = 1;
        ant_number = 2000;
        init_ant_variance = 2;
        periodic_boundaries = true;

        % +++++++++++++++++++++++++++++
        %             Ant model
        % +++++++++++++++++++++++++++++
        %------parameters-------
        ant_velocity = 1;
        nest_drift_velocity = .1;
        gradient_coupling = .1;
        beta = 1;
        noise_strength = 1;
        reversal_strength = 5;
        obstacle_strength = 5;
        %-----dynamical terms----
        gradient_force = true;
        nest_attraction_with_food = true;
        reversal_after_food = false;
        ballistic_coupling = true;
        reverse_at_food = true;
        angle_max = pi/20;
        gamma = 2;
        % +++++++++++++++++++++++++++++
        %         Pheromone model
        % +++++++++++++++++++++++++++++
        evap_rate = .001;
        deposition_rate = .02;
        max_pheromone = false;
        free_deposition_rate = 0;
        release_delay = 10;

        % +++++++++++++++++++++++++++++
        %             Geometry
        % +++++++++++++++++++++++++++++
        Lx = 100;
        Ly = 100;
        %----nest----
        nest_center =  [10 50];
        nest_radius = 4;
        %----food----
        food_center = [60 50];
        food_radius = 6;
        food_boundary_center = [0 0];
        food_boundary_radius = 10;
        %----obstacles-----
        obstacles = [];
        
        % +++++++++++++++++++++++++++++
        %         Visualisation
        % +++++++++++++++++++++++++++++
        show_figure = true;
        show_mesh = false;
        show_quiver = false;
        show_pheromones = false;
        display_interval = 10; % timesteps to skip between visualizations.
        display_pheromones_interval = 100; % timesteps to skip between pheromones visualizations.
        pause_time = .004;
        %-------movies-------
        create_movie = false;
        movie_num_frames = 200;
        movie_fps = 10;
        movie_filename = 'ant_sim';
        movie_quality = 100;
        
        % +++++++++++++++++++++++++++++
        %         Output
        % +++++++++++++++++++++++++++++
        
    end
    
    methods
        function obj = SimParams()
        end
    end
    
end

