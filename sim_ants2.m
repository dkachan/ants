function  out  = sim_ants2( p )
%SIM_ANTS Simulates ants with parameters p.
%   Detailed explanation goes here
% ********************************************************************
% *****                     Parameters                           *****
% ********************************************************************
%++++++++++++++++++++++++++++++
%       SIM SETTINGS
% +++++++++++++++++++++++++++++
N_time_steps = p.N_time_steps;
ant_number = p.ant_number;
dt = p.dt;

% +++++++++++++++++++++++++++++
%    GEOMETRY
% +++++++++++++++++++++++++++++
nest_center = p.nest_center;

% +++++++++++++++++++++++++++++
%    INITIAL ANT DISTRIBUTION
% +++++++++++++++++++++++++++++
init_ant_pos = ones(ant_number,1)*nest_center+(rand(ant_number,2)-0.5)*p.init_ant_variance;

% +++++++++++++++++++++++++++++
%            ANT MODEL
% +++++++++++++++++++++++++++++

ant_velocity = p.ant_velocity;
noise_strength = p.noise_strength;

% +++++++++++++++++++++++++++++
%    PREALLOCATION
% +++++++++++++++++++++++++++++

ant_pos = init_ant_pos;
prev_angle = 2*pi*rand(ant_number, 1);
delta_theta = zeros(ant_number,1);
trail_angle =p.trail_angle;
nsample = 0;
%--------output------------
%food_trace = zeros(N_time_steps,1);

% +++++++++++++++++++++++++++++
%    VISUALIZATION
% +++++++++++++++++++++++++++++
show_figure = p.show_figure;
pause_time = p.pause_time;
create_movie = p.create_movie;
display_interval = p.display_interval;

% ********************************************************************
% *****                      MAIN FLOW                           *****
% ********************************************************************


if create_movie
    frame_interval = floor(N_time_steps/p.movie_num_frames);
    vidObj = VideoWriter(p.movie_filename);
    vidObj.Quality = p.movie_quality;
    vidObj.FrameRate = p.movie_fps;
    open(vidObj);
    movie_counter=1;
end

if show_figure
    figure(1)
    clf(1)
    hold on
    xlim([0 1])
    ylim([0 1])
end


%-----Pre-compute pheromone map (static case)------
n_ptx = 100;
n_pty  =100;
[X, Y]=meshgrid(linspace(0,1,n_ptx),linspace(0,1,n_pty));

i = 1;
while (1)
    [x(i), y(i)] = ginput (1) ;
    if x(i)>0;
        plot (x(i), y(i), 'o');
    end
    if x(i)<0;
        x(i)=[];
        y(i)=[];
        break; end;
    i = i + 1;
end
curvexy(:,1)=x;
curvexy(:,2)=y;

mapxy = [X(:) Y(:)];
[xy,distance_to_trail_map_vec,~] = distance2curve(curvexy,mapxy,'linear');
trail_angle_map_vec = atan((Y(:)-xy(:,2))./(X(:)-xy(:,1)))+pi/2;

%-----Visualize if necessary---------
if show_figure
    figure(1)
    clf(1)
    quiver(X(:),Y(:),cos(trail_angle_map_vec)./max(0.01,distance_to_trail_map_vec),sin(trail_angle_map_vec)./max(0.01,distance_to_trail_map_vec))
    hold on
    figure(2)
    clf(2)
    hold on
end


%-----Time Iteration--------

for it = 1:N_time_steps
    % calculate the delta between currant angle and trail angle
    
    ix = floor(n_ptx*(ant_pos(:,1)-mod(ant_pos(:,1),1/n_ptx))+1e-3);
    iy = floor(n_pty*(ant_pos(:,2)-mod(ant_pos(:,2),1/n_pty))+1e-3);
    
    ix(ix==0) = 1;
    iy(iy==0) = 1;
    
    trail_angle = trail_angle_map_vec(100*(ix-1)+iy);
    distance_to_trail = distance_to_trail_map_vec(100*(ix-1)+iy);

    delta_angle = min(abs(mod(prev_angle-trail_angle,pi)),abs(mod(prev_angle-trail_angle,-pi)));
    ant_pheromones = 1./(max(delta_angle,0.01).*((distance_to_trail).^(0.75)));
   
    %Dynamics.
    %----angular noise----
    tanh_factor = tanh(p.beta*ant_pheromones);
    delta_theta = noise_strength*(1-tanh_factor).*(rand(ant_number, 1) - .5);
    %take 3
      
    if it < N_time_steps
        new_angle = prev_angle + delta_theta;
        ant_pos = ant_pos + ant_velocity*[cos(new_angle) sin(new_angle)]*dt;
        prev_angle = mod(new_angle,2*pi);
    end
    
    %------periodic_bc------
    if p.periodic_boundaries
        ant_pos((ant_pos(:,1)>1),1) = ant_pos((ant_pos(:,1)>1),1) - 1;
        ant_pos((ant_pos(:,1)<0),1) = ant_pos((ant_pos(:,1)<0),1) + 1;
        ant_pos((ant_pos(:,2)>1),2) = ant_pos((ant_pos(:,2)>1),2) - 1;
        ant_pos((ant_pos(:,2)<0),2) = ant_pos((ant_pos(:,2)<0),2) + 1;
    end
    %----------Real time visualization-----------
    if  show_figure && mod(it,display_interval) == 0
        figure(1)
        ant_loc_plot = plot(ant_pos(abs(delta_angle)>0,1),ant_pos(abs(delta_angle)>0,2),'ro');
        ant_pos_in = ant_pos(distance_to_trail<0.02,:);
        delta_angle_in = delta_angle(distance_to_trail<0.02);
        ant_loc_plot2 = plot(ant_pos_in(abs(delta_angle_in)<0.1,1),ant_pos_in(abs(delta_angle_in)<0.1,2),'g+');
        ant_loc_plot1 = plot(ant_pos(1,1),ant_pos(1,2),'m*','markersize',10);
        if create_movie && mod(it,frame_interval) == 0
           writeVideo(vidObj,getFrame(gca))
           movie_counter=movie_counter+1;
        end
        xlim([0,1]);
        ylim([0,1]);
        pause(pause_time)
        if it<N_time_steps
            delete(ant_loc_plot)
            delete(ant_loc_plot2)
            delete(ant_loc_plot1)
        end
    end
    
%     if mod(it,display_interval) == 0
%         nsample = nsample+1;
%         recruted_ants = length(ant_pos_in(abs(delta_angle_in)<0.1,1));
%         figure(2)
%         hold on
%         plot(nsample,recruted_ants,'r+')
%         xlim([0,nsample]);
%         ylim([0,recruted_ants+100]);
%     end
    

end


%--------Close video object-------------
if create_movie
    close(vidObj)
end

%placeholder
out = 1;

end