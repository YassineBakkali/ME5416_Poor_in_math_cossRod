clc;
clear all;
close all;

% Get the current working directory
currentFolder = pwd;
relativePath = fullfile(currentFolder, 'cossrod_funcs');
addpath(genpath(relativePath));
relativePath = fullfile(currentFolder, 'math_utils');
addpath(genpath(relativePath));
relativePath = fullfile(currentFolder, 'simulator');
addpath(genpath(relativePath));
relativePath = fullfile(currentFolder, 'actuators');
addpath(genpath(relativePath));
relativePath = fullfile(currentFolder, 'control');
addpath(genpath(relativePath));

L = 0.5; % Length of the rod [m]
r0 = 0.015; % radius of the rod [m]
rend = 0.0015;
start = [0; 0; 0];
nElements = 100; % Number of cross-sections
direction = [1.0; 0.0; 0.0];
normal = [0.0; 0.0; -1.0];

damping_cte = 10.0;

% Material Properties %

E0 = 1e5;               % Young modulus of rubber. [N/m^2]
nu = 0.5;              % Poisson ratio of rubber. [-]
G0 = E0 / (2*(1 + nu)); % Shear modulus of rubber. [N/m^2]

rho_material = 1000;    % Density of rubber        [kg/m^3]  
rho_medium = 1000;      % Density of water          [kg/m^3]

params = init_rod(nElements, direction, normal, L, r0, rend, rho_material, rho_medium, E0, start, G0, [], [], [], []);
dl = L / nElements;

dt = dl*0.05/2;
final_time = 0.5;
total_steps = final_time / dt;
coss_rod = cossrod_class(params);


% Creating the target %

cost_weight = struct('position', 0, 'director', 0);

target_cost_weight = struct('position', 1e6, 'director', 1e3);

target = [0.1; 0.15; 0.00];

dir = [0,0,1;0,1,0;-1,0,0];

targetPoint = TargetPoint(target,dir,cost_weight,target_cost_weight);

% Create muscles %

muscles = {};

LM_ratio_muscle_position = 0.0075/r0;

AN_ratio_radius = 0.002/r0;
TM_ratio_radius = 0.0045/r0;
LM_ratio_radius = 0.003/r0;

rod_area = pi .* coss_rod.radius.^2;
TM_rest_muscle_area = rod_area.*(TM_ratio_radius^2-AN_ratio_radius^2);
LM_rest_muscle_area = rod_area.*(LM_ratio_radius^2);

TM_max_muscle_stress =  15000.0;
LM_max_muscle_stress =  10000.0*20;

muscles{end+1} = transverse_muscle(nElements,TM_rest_muscle_area,TM_max_muscle_stress);
muscles{end+1} = longitudinal_muscle(nElements, pi*0.5*0, LM_ratio_muscle_position, LM_rest_muscle_area, LM_max_muscle_stress);
muscles{end+1} = longitudinal_muscle(nElements, pi*0.5*1, LM_ratio_muscle_position, LM_rest_muscle_area, LM_max_muscle_stress);
muscles{end+1} = longitudinal_muscle(nElements, pi*0.5*2, LM_ratio_muscle_position, LM_rest_muscle_area, LM_max_muscle_stress);
muscles{end+1} = longitudinal_muscle(nElements, pi*0.5*3, LM_ratio_muscle_position, LM_rest_muscle_area, LM_max_muscle_stress);

for idx = 1:numel(muscles)
    muscles{idx} = muscles{idx}.set_length_as_rest(coss_rod);
end

params = struct('stepsize', 1e-8, 'activation_diff_tol', 1e-10);

fwbw_algo = ForwardBackwardAlgo(coss_rod, params, muscles, targetPoint);

% Initialize the director as an identity matrix
director = eye(3);

% Calculate the vectors from the base and tip of the rod to the target
base_to_target = fwbw_algo.target.pos - coss_rod.pos_vects(:, 1);
tip_to_target = fwbw_algo.target.pos - coss_rod.pos_vects(:, end);

% Normalize the base_to_target and tip_to_target vectors
base_to_target = base_to_target / norm(base_to_target);
tip_to_target = tip_to_target / norm(tip_to_target);

% Calculate the new director vectors using cross products
director(2, :) = cross(base_to_target, tip_to_target);
director(1, :) = cross(director(2, :), tip_to_target);
director(3, :) = cross(director(1, :), director(2, :));

% Assign the computed director matrix to algo.objects and target.director_collection
fwbw_algo.target.dir = director;
targetPoint.dir(:, :, 1) = director;

fwbw_algo = fwbw_algo.runAlgo(100000);

% 
% simu = simulator_class(coss_rod, damping_cte, dt);
% simu = simu.integrate(final_time,total_steps);
% Assume simu.rod.pos_vects is your [3, n+1] matrix where each column is [x; y; z]
% and simu.rod.radius is a [1, n] vector containing the radius at each segment
x = fwbw_algo.rod.pos_vects(1, :);  % X coordinates
y = fwbw_algo.rod.pos_vects(2, :);  % Y coordinates
z = fwbw_algo.rod.pos_vects(3, :);  % Z coordinates
r = fwbw_algo.rod.radius;           % Radii for each segment

% Create a new figure
figure;

% Number of points around each circle
numPtsCirc = 20;

% Precompute sine and cosine for circle
theta = linspace(0, 2*pi, numPtsCirc);
sinTheta = sin(theta);
cosTheta = cos(theta);

% Hold on to add multiple plots
hold on;

% Plot each segment with its radius
for i = 1:length(x)-1
    % Direction from this point to the next
    direction = [x(i+1) - x(i); y(i+1) - y(i); z(i+1) - z(i)];
    direction = direction / norm(direction);  % Normalize

    % Find a normal vector using cross product
    normal = cross(direction, [1; 0; 0]); % Use a fixed vector for consistency
    normal = normal / norm(normal);       % Normalize

    % Another normal perpendicular to the first normal and direction
    binormal = cross(direction, normal);
    binormal = binormal / norm(binormal); % Normalize

    % Generate cylinder points
    circleX = x(i) + r(i) * cosTheta * normal(1) + r(i) * sinTheta * binormal(1);
    circleY = y(i) + r(i) * cosTheta * normal(2) + r(i) * sinTheta * binormal(2);
    circleZ = z(i) + r(i) * cosTheta * normal(3) + r(i) * sinTheta * binormal(3);

    % Generate the upper part of the cylinder
    upperX = x(i+1) + r(i) * cosTheta * normal(1) + r(i) * sinTheta * binormal(1);
    upperY = y(i+1) + r(i) * cosTheta * normal(2) + r(i) * sinTheta * binormal(2);
    upperZ = z(i+1) + r(i) * cosTheta * normal(3) + r(i) * sinTheta * binormal(3);

    % Use surf to plot the cylinder
    surf([circleX; upperX], ...
         [circleY; upperY], ...
         [circleZ; upperZ], ...
         'EdgeColor', 'none', 'FaceColor', 'interp');
end

% Set the axes limits to ensure all data is visible
xlim([min(x)-max(r), max(x)+max(r)]);
ylim([min(y)-max(r), max(y)+max(r)]);
zlim([min(z)-max(r), max(z)+max(r)]);

% Add labels and title for clarity
xlabel('X Axis');
ylabel('Y Axis');
zlabel('Z Axis');
title('3D Visualization of the Rod with Radii');

% Add grid for better visualization
grid on;

% Set the aspect ratio to equal to ensure all axes are scaled equally
axis equal;

% Set a 3D view angle
view(45, 25);

% End the hold on
hold off;